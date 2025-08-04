#!/usr/bin/env python3


from __future__ import annotations

import argparse, csv, random, sqlite3, sys, textwrap
from collections import Counter
from datetime import datetime
from multiprocessing import cpu_count, get_context
from pathlib import Path
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd

# ───────────────────────── Global tweaks ──────────────────────────

csv.field_size_limit(256 * 1024 * 1024)  # 256 MB field limit

# ───────────────────────── Helpers ────────────────────────────────

def log(msg: str):
    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}")

def comp(base: str) -> str:
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(base, base)

from numpy.random import SeedSequence, default_rng, Generator

def rng_streams(seed: Optional[int], n: int) -> List[Generator]:
    if seed is None:
        return [default_rng() for _ in range(n)]
    return [default_rng(s) for s in SeedSequence(seed).spawn(n)]

# ───────────────────────── CLI ────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--output")
    p.add_argument("--fasta_lengths_info")
    p.add_argument("--target_exome_region")
    p.add_argument("--set_seed", type=int)
    p.add_argument("--cores", type=int)
    p.add_argument("--mutational_context")
    p.add_argument("--Ts_Tv_ratio")
    p.add_argument("--variant_probability_table")
    p.add_argument("--sbs_signature_table")
    g1 = p.add_mutually_exclusive_group()
    g1.add_argument("--nr_copies", type=int)
    g1.add_argument("--multifasta_copies")
    g2 = p.add_mutually_exclusive_group()
    g2.add_argument("--mut_genome_fraction", type=float)
    g2.add_argument("--specific_mut_rate", type=float)
    p.add_argument("--fraction_positions", type=float)
    return p

# ──────────────────────── SQLite helpers ──────────────────────────

def sql_conn(header: str):
    return sqlite3.connect(Path("/usr/src/app/pipeline/SQL_database") / f"{header}.sqlite")

def fetch_rows(header: str, ids: Iterable[int]) -> pd.DataFrame:
    ids = sorted({int(i) for i in ids})
    if not ids:
        return pd.DataFrame()
    with sql_conn(header) as conn:
        ph = ",".join(["?"] * len(ids))
        df = pd.read_sql_query(f"SELECT * FROM {header} WHERE ROWID IN ({ph})", conn, params=ids)
    df["position"] = ids
    return df

# ───────── Variant‑probability table & sampler ────────────────────

def default_table(ts_prob: float, tv_prob: float) -> pd.DataFrame:
    rows = [(b, v, ts_prob if (b+v) in ("AG","GA","CT","TC") else tv_prob/2)
            for b in "ATCG" for v in "ATCG" if v != b]
    return pd.DataFrame(rows, columns=["Nucleotide","Variant","Probability"])

def build_sampler(tbl: pd.DataFrame):
    lut = {b:(g["Variant"].values,(g["Probability"].astype(float)/g["Probability"].sum()).values)
           for b,g in tbl.groupby("Nucleotide")}
    _rng=default_rng()
    def _s(base:str):
        if base=="N" or pd.isna(base):
            return base
        alleles,probs=lut[base]
        return _rng.choice(alleles,p=probs)
    return _s

# ───────── Context queries (pickle‑safe) ─────────────────────────

def ctx_query(header:str,col:str,rev:str,ctx:List[str],start:int,end:int):
    ph=",".join(["?"]*len(ctx))
    params=ctx+[start,end]+ctx+[start,end]
    sql=textwrap.dedent(f"""SELECT ROWID,'fwd' FROM {header} WHERE {col} IN ({ph}) AND ROWID BETWEEN ? AND ?
UNION ALL SELECT ROWID,'rev' FROM {header} WHERE {rev} IN ({ph}) AND ROWID BETWEEN ? AND ?""")
    with sql_conn(header) as conn:
        rows=conn.execute(sql,params).fetchall()
    fwd=[r[0] for r in rows if r[1]=='fwd']
    rev=[r[0] for r in rows if r[1]=='rev']
    return np.array(fwd,int),np.array(rev,int)

def ctx_batch(args):
    batch,hd,c1,c2,ctx=args
    return ctx_query(hd,c1,c2,ctx,int(batch.iloc[0,1]),int(batch.iloc[-1,2]))

# ───────── SBS helpers ───────────────────────────────────────────

def load_sbs(path:Path):
    df=pd.read_csv(path)
    parts=df["X"].str.extract(r"^(?P<a>[A-Z]).*\[(?P<b>[A-Z])>(?P<v>[A-Z])\].*(?P<c>[A-Z])")
    df["context"]=parts["a"]+parts["b"]+parts["c"]
    df["variant"]=parts["v"]
    tf=df[["context","variant",df.columns[1]]]
    tf.columns=["context","variant","occurance_prob"]
    ctx_sum=tf.groupby("context",as_index=False)["occurance_prob"].sum()
    pivot=(tf.pivot_table(index="context",columns="variant",values="occurance_prob",aggfunc="sum").fillna(0.0).reset_index())
    return ctx_sum,pivot

def sbs_count(args):
    batch,hd,contexts=args
    ph=",".join(["?"]*len(contexts))
    loc=Counter()
    with sql_conn(hd) as conn:
        for _,(_,s,e) in batch.iterrows():
            q=textwrap.dedent(f"SELECT tri_context,tri_context_rev_comp,COUNT(*) FROM {hd} WHERE (tri_context IN ({ph}) OR tri_context_rev_comp IN ({ph})) AND ROWID BETWEEN ? AND ? GROUP BY tri_context,tri_context_rev_comp")
            params=contexts+contexts+[int(s),int(e)]
            for c1,c2,cnt in conn.execute(q,params):
                loc[c1]+=cnt; loc[c2]+=cnt
    return loc

# ───────── CSV utilities ─────────────────────────────────────────

def wcsv(df:pd.DataFrame,path:Path):
    df.to_csv(path,sep=";",index=False,quoting=csv.QUOTE_NONE,lineterminator="\n")

def cat(parts:List[Path],dst:Path):
    frames=[pd.read_csv(p,sep=";") for p in parts]
    pd.concat(frames,ignore_index=True).to_csv(dst,sep=";",index=False,quoting=csv.QUOTE_NONE,lineterminator="\n")

def add_suffix(src:Path,dst:Path):
    df=pd.read_csv(src,sep=";")
    df.iloc[:,0]=[f"{n}_MG{i+1}" for i,n in enumerate(df.iloc[:,0])]
    wcsv(df,dst)

# ───────── Misc helpers ─────────────────────────────────────────

def split_df(df:pd.DataFrame,n:int):
    if n<=1 or df.empty:
        return [df]
    k,m=divmod(len(df),n)
    return [df.iloc[i*k+min(i,m):(i+1)*k+min(i+1,m)] for i in range(n)]

def mg_row(idx:int,posvar:np.ndarray,hd:str,rng:Generator):
    if posvar.size==0:
        return {"genome_name":f"{hd}_MG{idx+1}","variant_position":"No_mutations","copy_number":1}
    n=rng.integers(1,posvar.size+1)
    sel=rng.choice(posvar,size=n,replace=False)
    sel=sorted(sel,key=lambda s:int("".join(filter(str.isdigit,s))))
    return {"genome_name":f"{hd}_MG{idx+1}","variant_position":",".join(sel),"copy_number":1}

# ───────── Core per‑FASTA processor ────────────────────────────

def run_fasta(row:pd.Series,bed:pd.DataFrame,a,samp,rng,out:Path,cores:int,sbs):
    hd=row["fasta_names"]
    log(f"→ Processing {hd}")
    log(f"→ {hd}")
    bed_chr=bed[bed.iloc[:,0]==hd]
    if bed_chr.empty:
        log("(no BED rows)"); return
    ctx=a.mutational_context.split(",") if a.mutational_context else None
    if ctx:
        L=len(ctx[0]); col,rev={1:("nucleotide","comp_nucleotide"),3:("tri_context","tri_context_rev_comp"),5:("penta_context","penta_context_rev_comp")}[L]
    else:
        col,rev="nucleotide","comp_nucleotide"
    batches=split_df(bed_chr,cores)
    if ctx:
        with get_context("fork").Pool(cores) as p:
            res=p.map(ctx_batch,[(b,hd,col,rev,ctx) for b in batches])
        pos5=np.unique(np.concatenate([r[0] for r in res])); pos3=np.unique(np.concatenate([r[1] for r in res]))
        all_pos=np.unique(np.concatenate([pos5,pos3]))
    else:
        all_pos=np.hstack([np.arange(int(s),int(e)+1) for _,(_,s,e) in bed_chr.iterrows()])
        pos5=pos3=np.array([],int)

    # Mode A
    if a.mut_genome_fraction is not None and a.specific_mut_rate is None and sbs is None:
        log(f"→ Mode A selected: mut_genome_fraction={a.mut_genome_fraction}")
        pos_fraction=a.fraction_positions or 0.05
        n_mut=int(round(len(all_pos)*pos_fraction))
        picked=rng.choice(all_pos,size=n_mut,replace=False); picked.sort()
        if ctx:
            pick5=picked[np.isin(picked,pos5)]; pick3=picked[np.isin(picked,pos3)]
        else:
            half=rng.choice(picked,size=len(picked)//2,replace=False); pick5,pick3=half,np.setdiff1d(picked,half)
        mut_tbl=fetch_rows(hd,picked)
        mut_tbl["Variant"]=mut_tbl["nucleotide"].apply(samp)
        mut_tbl.loc[mut_tbl["position"].isin(pick3),"Variant"]=mut_tbl.loc[mut_tbl["position"].isin(pick3),"Variant"].map(comp)
        mut_tbl.sort_values("position",inplace=True)
        wcsv(mut_tbl,out/f"{hd}_mutations_table.csv")

        copies=row["copy_number"]; MGs=int(round(copies*a.mut_genome_fraction))
        NMG=pd.DataFrame({"genome_name":[f"{hd}_NMG"],"variant_position":["No_mutations"],"copy_number":[copies-MGs]})
        NMG_tmp=out/f"{hd}_NMG_tmp.csv"; wcsv(NMG,NMG_tmp)
        posvar=(mut_tbl["position"].astype(str)+mut_tbl["Variant"].astype(str)).values
        streams=rng_streams(a.set_seed,MGs)
        with get_context("fork").Pool(cores) as p:
            mg_rows=p.starmap(mg_row,[(i,posvar,hd,streams[i]) for i in range(MGs)])
        MG_tmp=out/f"{hd}_MG_tmp.csv"; wcsv(pd.DataFrame(mg_rows),MG_tmp)
        MG_tmp2=out/f"{hd}_MG_tmp_renamed.csv"; add_suffix(MG_tmp,MG_tmp2)
        cat([NMG_tmp,MG_tmp2], out/f"{hd}_sample_table.csv")
        NMG_tmp.unlink(); MG_tmp.unlink(); MG_tmp2.unlink();
        return

    # Mode B
    if a.specific_mut_rate is not None:
        log(f"→ Mode B selected: specific_mut_rate={a.specific_mut_rate}")
        ex_len=(bed_chr.iloc[:,2]-bed_chr.iloc[:,1]+1).sum()
        lam=a.specific_mut_rate*ex_len
        copies=row["copy_number"]
        pois=rng.poisson(lam,size=copies); has=pois>0
        if not has.any():
            wcsv(pd.DataFrame({"genome_name":[f"{hd}_NMG"],"variant_position":["No_mutations"],"copy_number":[copies]}), out/f"{hd}_sample_table.csv"); return
        n_mut=pois[has]; MGn=len(n_mut)
        MG_df=pd.DataFrame({"genome_name":[f"{hd}_MG{i+1}" for i in range(MGn)],"nr_mutations":n_mut})
        pick=lambda k:rng.choice(all_pos,size=k,replace=False).tolist()
        MG_df["positions"]=MG_df["nr_mutations"].apply(pick)
        picked_all=np.unique(np.concatenate(MG_df["positions"].values))
        if ctx:
            pick5=picked_all[np.isin(picked_all,pos5)]; pick3=picked_all[np.isin(picked_all,pos3)]
        else:
            half=rng.choice(picked_all,size=len(picked_all)//2,replace=False); pick5,pick3=half,np.setdiff1d(picked_all,half)
        mut_tbl=fetch_rows(hd,picked_all)
        mut_tbl["Variant"]=mut_tbl["nucleotide"].apply(samp)
        mut_tbl.loc[mut_tbl["position"].isin(pick3),"Variant"]=mut_tbl.loc[mut_tbl["position"].isin(pick3),"Variant"].map(comp)
        mut_tbl.sort_values("position",inplace=True); wcsv(mut_tbl,out/f"{hd}_mutations_table.csv")
        pv_map={p:v for p,v in zip(mut_tbl["position"],mut_tbl["Variant"])}
        MG_df["variant_position"]=MG_df["positions"].apply(lambda lst:",".join(f"{p}{pv_map[p]}" for p in sorted(lst)))
        sample_MG=MG_df.groupby("variant_position",as_index=False).agg(genome_name=("genome_name","first"), copy_number=("variant_position","size"))
        NMG=pd.DataFrame({"genome_name":[f"{hd}_NMG"],"variant_position":["No_mutations"],"copy_number":[copies-MGn]})
        wcsv(pd.concat([NMG,sample_MG],ignore_index=True), out/f"{hd}_sample_table.csv")
        return

    # Mode C
    if sbs is not None and a.mut_genome_fraction is not None:
        log("→ Mode C selected: SBS signature mode")
        ctx_sum,pivot=sbs; contexts=ctx_sum["context"].tolist()
        with get_context("fork").Pool(cores) as p:
            counts=p.map(sbs_count,[(b,hd,contexts) for b in batches])
        occ=Counter(); [occ.update(c) for c in counts]
        ctx_df=pd.DataFrame({"context":list(occ.keys()),"genome_occ":list(occ.values())})
        ctx_df=ctx_df.merge(ctx_sum,on="context",how="left").fillna({"occurance_prob":0})
        ctx_df["prob"]=ctx_df["genome_occ"]/ctx_df["genome_occ"].sum()
        ctx_df["adj"]=ctx_df["occurance_prob"]*ctx_df["prob"]
        total_occ=ctx_df["genome_occ"].sum(); pos_fraction=a.fraction_positions or 0.05
        n_mut=int(round(total_occ*pos_fraction))
        draws=rng.multinomial(n_mut, ctx_df["adj"]/ctx_df["adj"].sum()); ctx_df["draw"]=draws
        def sample_pos(row):
            if row.draw==0: return []
            acc=[]
            with sql_conn(hd) as conn:
                for _,(_,s,e) in bed_chr.iterrows():
                    q=f"SELECT ROWID FROM {hd} WHERE (tri_context=? OR tri_context_rev_comp=?) AND ROWID BETWEEN ? AND ?"
                    cur=conn.execute(q,(row.context,row.context,int(s),int(e)))
                    acc.extend([r[0] for r in cur.fetchall()])
            if not acc: return []
            return rng.choice(acc,size=row.draw,replace=False).tolist()
        ctx_df["pos"]=ctx_df.apply(sample_pos,axis=1)
        sampled=np.unique(np.concatenate(ctx_df["pos"].values))
        mut_tbl=fetch_rows(hd,sampled)
        prob_map={r.context:r[["A","C","G","T"]].values.astype(float) for _,r in pivot.iterrows()}
        bases=np.array(list("ACGT"))
        def assign(row):
            tri,tri_rev=row.tri_context,row.tri_context_rev_comp
            if tri in prob_map:
                p=prob_map[tri]/prob_map[tri].sum(); return rng.choice(bases,p=p)
            if tri_rev in prob_map:
                p=prob_map[tri_rev]/prob_map[tri_rev].sum(); return comp(rng.choice(bases,p=p))
            return np.nan
        mut_tbl["Variant"]=mut_tbl.apply(assign,axis=1)
        mut_tbl.sort_values("position",inplace=True); wcsv(mut_tbl,out/f"{hd}_mutations_table.csv")
        copies=row["copy_number"]; MGs=int(round(copies*a.mut_genome_fraction))
        NMG=pd.DataFrame({"genome_name":[f"{hd}_NMG"],"variant_position":["No_mutations"],"copy_number":[copies-MGs]})
        NMG_tmp=out/f"{hd}_NMG_tmp.csv"; wcsv(NMG,NMG_tmp)
        posvar=(mut_tbl["position"].astype(str)+mut_tbl["Variant"].astype(str)).values
        streams=rng_streams(a.set_seed,MGs)
        with get_context("fork").Pool(cores) as p:
            mg_rows=p.starmap(mg_row,[(i,posvar,hd,streams[i]) for i in range(MGs)])
        MG_tmp=out/f"{hd}_MG_tmp.csv"; wcsv(pd.DataFrame(mg_rows),MG_tmp)
        MG_tmp2=out/f"{hd}_MG_tmp_renamed.csv"; add_suffix(MG_tmp,MG_tmp2)
        cat([NMG_tmp,MG_tmp2], out/f"{hd}_sample_table.csv")
        NMG_tmp.unlink(); MG_tmp.unlink(); MG_tmp2.unlink();
        return

    sys.exit("Bad flag combination – check parameters")

# ───────── Main entry ───────────────────────────────────────────

def main():
    a=build_parser().parse_args()
    log("============================================")
    log("  GENOMICON-Seq Mutation Simulation Started")
    log("============================================")
    log(f"→ Output directory: {a.output}")
    log(f"→ FASTA info file: {a.fasta_lengths_info}")
    log(f"→ BED region file: {a.target_exome_region}")
    log(f"→ Seed: {a.set_seed if a.set_seed is not None else 'None'}")
    log(f"→ Requested cores: {a.cores if a.cores is not None else 'None'}")
    log(f"→ Mutational context: {a.mutational_context if a.mutational_context else 'None'}")
    log(f"→ Transition/Transversion ratio: {a.Ts_Tv_ratio if a.Ts_Tv_ratio else 'None'}")
    log(f"→ Variant probability table: {a.variant_probability_table if a.variant_probability_table else 'None'}")
    log(f"→ SBS signature table: {a.sbs_signature_table if a.sbs_signature_table else 'None'}")
    log(f"→ Number of genome copies: {a.nr_copies if a.nr_copies else ('[multifasta]' if a.multifasta_copies else 'Default 100')}")
    log(f"→ Mutation genome fraction: {a.mut_genome_fraction if a.mut_genome_fraction is not None else 'None'}")
    log(f"→ Specific mutation rate: {a.specific_mut_rate if a.specific_mut_rate is not None else 'None'}")
    log(f"→ Fraction of positions to mutate: {a.fraction_positions if a.fraction_positions is not None else 'Default 0.05'}")

   

    for req in ("output","fasta_lengths_info","target_exome_region"):
        if getattr(a,req) is None:
            sys.exit(f"--{req} was not provided")
    out=Path(a.output).resolve(); out.mkdir(parents=True,exist_ok=True)
    if a.set_seed is not None and a.set_seed<0:
        sys.exit("Seed must be ≥0")
    if a.set_seed is not None:
        random.seed(a.set_seed); np.random.seed(a.set_seed)
    rng = np.random.default_rng(a.set_seed)

    avail = cpu_count()
    requested = a.cores if a.cores is not None else 1

    if requested > avail:
        fallback = max(avail - 2, 1)
        log(f"Warning: requested {requested} cores, but only {avail} available. Using {fallback} instead.")
        cores = fallback
    else:
        cores = requested
    log(f"Cores requested: {requested} | Cores available: {avail} | Cores used: {cores}")

    fasta_df = pd.read_csv(a.fasta_lengths_info)
    bed_df=pd.read_csv(a.target_exome_region,sep="\t",header=None)
    log("Input tables loaded")
    if a.nr_copies:
        if a.nr_copies<=0: sys.exit("--nr_copies must be positive")
        fasta_df["copy_number"]=a.nr_copies
    elif a.multifasta_copies:
        fasta_df = fasta_df.merge(
            pd.read_csv(a.multifasta_copies).rename(columns={"fasta_name": "fasta_names"}),
            on="fasta_names"
        )
    else:
        log("Copy number not provided ⇒ default 100"); fasta_df["copy_number"]=100
    if a.Ts_Tv_ratio:
        if ":" not in a.Ts_Tv_ratio: sys.exit("--Ts_Tv_ratio must be like 3:1")
        ts,tv=map(int,a.Ts_Tv_ratio.split(":")); tot=ts+tv
        var_tbl=default_table(ts/tot,tv/tot)
    elif a.variant_probability_table:
        path=Path(a.variant_probability_table).resolve()
        if not str(path).startswith("/usr/src/app/pipeline"):
            sys.exit("variant_probability_table must sit under /usr/src/app/pipeline")
        var_tbl=pd.read_csv(path)
    else:
        var_tbl=default_table(0.5,0.5)
    sampler=build_sampler(var_tbl)
    if a.sbs_signature_table:
        sbs_tables=load_sbs(Path(a.sbs_signature_table))
    else:
        sbs_tables=None
    for _,r in fasta_df.iterrows():
        run_fasta(r,bed_df,a,sampler,rng,out,cores,sbs_tables)
    
    log("============================================")
    log("Mutation simulation complete for all entries.")
    log("============================================")



if __name__=="__main__":
    main()
