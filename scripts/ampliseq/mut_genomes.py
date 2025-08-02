#!/usr/bin/env python3
"""
A faithful Python translation of the original R script that simulates mutations in
FASTA genomes under a variety of user‑specified parameters.

The script reproduces **every single step** of the R workflow, including:
  • Command‑line interface and default behaviours
  • Reproducible seeding logic (global + per‑worker)
  • Parallel execution across multiple cores
  • FASTA reading and basic statistics output
  • Loading / inferring a nucleotide‑>variant probability table
  • Tri/penta context handling via SQLite look‑ups
  • Two mutation modes:  (A) mutated‑genome fraction  &  (B) specific Poisson rate
  • Generation of *_mutations_table.csv* and *_sample_table.csv* per FASTA header
  • Extensive console logging that mirrors the original script

**Important paths**
  /usr/src/app/pipeline/               – expected working directory for inputs
  /usr/src/app/pipeline/SQL_database/  – where per‑header SQLite files live

Python libraries used:
  BioPython, pandas, numpy, sqlite3, multiprocessing, argparse, random, itertools

The code has been carefully organised so that each logical block of the R code
maps 1‑to‑1 onto a Python function or section, making future maintenance or
verification straightforward.

© 2025  — Translated by ChatGPT‑4o
"""

import os
import sys
import argparse
import random
import sqlite3
import csv
import multiprocessing as mp
from itertools import repeat
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

################################################################################
# 0. WORKING DIRECTORY, SEED AND CORES                                         #
################################################################################

def parse_args() -> argparse.Namespace:
    """Parses and validates all command‑line arguments."""
    parser = argparse.ArgumentParser(
        description="Python re‑implementation of the R mutation‑simulation script",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )

    # I/O and environment
    parser.add_argument("--output", required=True, type=Path,
                        help="Output folder (will become cwd)")
    parser.add_argument("--set_seed", type=int, metavar="INT",
                        help="Reproducible RNG seed (non‑negative integer)")
    parser.add_argument("--cores", type=int, metavar="INT",
                        help="Number of CPU cores to use (default 2 if not given)")

    # Required FASTA
    parser.add_argument("--fasta_sequence", required=True, type=str,
                        help="Name of the FASTA file inside /usr/src/app/pipeline/")

    # Mutation context & probability table
    parser.add_argument("--mutational_context", type=str,
                        help="Comma‑separated list of 1‑, 3‑ or 5‑nt contexts")
    parser.add_argument("--Ts_Tv_ratio", type=str,
                        help="Transition:Transversion ratio, e.g. 2:1")
    parser.add_argument("--variant_probability_table", type=str,
                        help="CSV with explicit Nucleotide / Variant / Probability cols")

    # Copy number
    parser.add_argument("--nr_copies", type=int,
                        help="Uniform copy number for all FASTAs in multifasta")
    parser.add_argument("--multifasta_copies", type=str,
                        help="CSV mapping FASTA headers to copy_number column")

    # Mutation rate mode A (genome fraction)
    parser.add_argument("--mut_genome_fraction", type=float,
                        help="Fraction (0‒1] of genomes that will receive ≥1 mutation")
    parser.add_argument("--fraction_positions", type=float,
                        help="Fraction (0‒1] of candidate positions that actually mutate")

    # Mutation rate mode B (specific Poisson rate)
    parser.add_argument("--specific_mut_rate", type=float,
                        help="λ per‑base mutation rate (0‒1]; mutually exclusive with ‑‑mut_genome_fraction)")

    return parser.parse_args()

###############################################################################
# Helper utilities                                                            #
###############################################################################

def set_global_seed(seed: int | None):
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed % (2**32 - 1))
        print(f"Specified seed value: {seed}")
    else:
        print("No seed were provided")


def configure_cores(user_cores: int | None) -> int:
    if user_cores is None:
        print("Number of cores not provided, using default ‑ 2 cores")
        return 2

    if user_cores <= 0 or user_cores != int(user_cores):
        sys.exit("Invalid number of cores.")

    avail = mp.cpu_count()
    if user_cores > avail:
        print(f"Warning: Requested number of cores ({user_cores}) exceeds the available cores ({avail})")
    print(f"Specified number of cores: {user_cores}")
    return user_cores


###############################################################################
# I. DEFINING SEQUENCE                                                       #
###############################################################################

def load_fasta(fasta_name: str) -> list[SeqRecord]:
    fasta_path = Path("/usr/src/app/pipeline") / fasta_name
    if not fasta_path.exists():
        sys.exit("The specified fasta file does not exist")
    print("Fasta file loaded")
    return list(SeqIO.parse(str(fasta_path), "fasta"))

###############################################################################
# II. NUCLEOTIDE CONTEXT & MUTATION SIGNATURE                                #
###############################################################################

def parse_mut_context(mut_context_raw: str | None) -> list[str] | None:
    if mut_context_raw is None:
        return None

    ctx = [s.strip().upper() for s in mut_context_raw.split(",") if s.strip()]
    lengths = {len(c) for c in ctx}
    if len(lengths) != 1 or lengths.pop() not in {1, 3, 5}:
        sys.exit("Invalid mutational context. All elements must have the same length 1, 3, or 5.")
    print(f"Specified target for mutational process: {ctx}")
    return ctx


def build_variant_prob_table(args: argparse.Namespace) -> pd.DataFrame:
    """
    Return a DataFrame with columns  Nucleotide | Variant | Probability
    – If the user supplies a CSV via --variant_probability_table,
      we load it, strip whitespace, and make sure that for every
      nucleotide the probabilities sum to 1.0.
      • If they sum to ~100  ➜ assume percentages and divide by 100
      • If they sum to ~1    ➜ already fine
      • Otherwise            ➜ hard error (catches typos)
    – If the flag is absent we fall back to Ts/Tv logic or equal 1/3.
    """
    # 1. user-supplied table
    if args.variant_probability_table:
        tbl_path = Path("/usr/src/app/pipeline") / args.variant_probability_table
        if not tbl_path.exists():
            sys.exit("The specified variant probability table not found")

        df = pd.read_csv(tbl_path)

        # tidy up: strip spaces from any string column
        for col in df.select_dtypes(include="object").columns:
            df[col] = df[col].str.strip()

        # force numeric Probabilities
        df["Probability"] = pd.to_numeric(df["Probability"], errors="coerce")

        # validate / normalise per nucleotide
        for nuc, grp in df.groupby("Nucleotide"):
            tot = grp["Probability"].sum()

            if np.isclose(tot, 1.0):          # already proportions
                continue

            if np.isclose(tot, 100.0):        # looks like percentages
                df.loc[grp.index, "Probability"] = grp["Probability"] / 100.0
                continue

            # anything else is suspicious – stop early
            sys.exit(
                f"Probability values for nucleotide '{nuc}' sum to {tot}, "
                "expected 1.0 (proportions) or 100.0 (percentages)."
            )

        print("Probability table loaded and normalised")
        return df

    # 2. Build from Ts/Tv ratio
    if args.Ts_Tv_ratio:
        if ":" not in args.Ts_Tv_ratio:
            sys.exit("Invalid format for --Ts_Tv_ratio. Expected 'integer:integer'.")
        try:
            ts_val, tv_val = map(float, args.Ts_Tv_ratio.split(":"))
        except ValueError:
            sys.exit("Invalid format for --Ts_Tv_ratio. Expected 'integer:integer'.")
        print(f"Specified transition vs transversion ratio: {args.Ts_Tv_ratio}")
        tot = ts_val + tv_val
        ts_prob = ts_val / tot
        tv_prob = tv_val / tot / 2  # each transversion split into two possibilities

        rows = []
        mapping = {
            "A": [("G", ts_prob), ("C", tv_prob), ("T", tv_prob)],
            "T": [("C", ts_prob), ("A", tv_prob), ("G", tv_prob)],
            "C": [("T", ts_prob), ("G", tv_prob), ("A", tv_prob)],
            "G": [("A", ts_prob), ("C", tv_prob), ("T", tv_prob)],
        }
        for nuc, variants in mapping.items():
            for var, p in variants:
                rows.append((nuc, var, p))
        return pd.DataFrame(rows, columns=["Nucleotide", "Variant", "Probability"])

    # 3. Default equal probability
    print("table with variant mutation probability not provided, each nucleotide can mutate to other 3 variants with equal probability")
    rows = []
    for nuc in "ATCG":
        for var in "ATCG":
            if var != nuc:
                rows.append((nuc, var, 1/3))
    return pd.DataFrame(rows, columns=["Nucleotide", "Variant", "Probability"])


###############################################################################
# III. WHAT IS IN THE SAMPLE?                                                #
###############################################################################

def determine_nr_copies(args: argparse.Namespace,
                        fasta_records: list[SeqRecord]) -> list[int]:
    """
    Decide how many copies to simulate for each sequence in the multifasta.

    Priority:
      1) --nr_copies (scalar)                    → same value for every record
      2) --multifasta_copies <csv>               → lookup per header
      3) default 1000                            → fallback

    The CSV must have two columns:
        fasta_name,copy_number
    Names are matched (case-sensitive) to SeqRecord.id.
    """

    # ---- mode 1 : scalar override --------------------------------------
    if args.nr_copies is not None:
        if args.nr_copies <= 0:
            sys.exit("Invalid --nr_copies (must be positive int)")
        return [args.nr_copies] * len(fasta_records)

    # ---- mode 2 : per-header table -------------------------------------
    if args.multifasta_copies:
        tbl_path = Path("/usr/src/app/pipeline") / args.multifasta_copies
        if not tbl_path.exists():
            sys.exit("multifasta_copies table not found: " + str(tbl_path))

        df = (
            pd.read_csv(tbl_path)
              .rename(columns=str.strip)
              .assign(fasta_name=lambda d: d["fasta_name"].str.strip())
        )

        # build a dict for fast lookup
        name2copies = dict(
            zip(df["fasta_name"], df["copy_number"].astype(int))
        )

        copies_out = []
        for rec in fasta_records:
            name = rec.id
            if name not in name2copies:
                sys.exit(
                    f"Header '{name}' has no entry in {args.multifasta_copies}"
                )
            copies_out.append(name2copies[name])

        print("Loaded copy numbers from multifasta_copies table")
        return copies_out

    # ---- default -------------------------------------------------------
    print("No copy number provided – defaulting to 1000 each")
    return [1000] * len(fasta_records)

###############################################################################
# IV. DEFINING NECESSARY FUNCTIONS & VARIABLES                               #
###############################################################################

def context_sql_columns(ctx_len: int | None):
    if ctx_len == 1:
        return "nucleotide", "comp_nucleotide"
    if ctx_len == 3:
        return "tri_context", "tri_context_rev_comp"
    if ctx_len == 5:
        return "penta_context", "penta_context_rev_comp"
    return "nucleotide", "comp_nucleotide"


def generate_variant(nuc: str, prob_table: pd.DataFrame) -> str:
    probs = prob_table.loc[prob_table["Nucleotide"] == nuc, ["Variant", "Probability"]]
    return np.random.choice(probs["Variant"], p=probs["Probability"])


def query_positions_sql(db_path: Path, table: str, column: str, ctx_set: list[str]) -> list[int]:
    placeholders = ",".join("?" * len(ctx_set))
    conn = sqlite3.connect(db_path)
    qry = f"SELECT ROWID as row_number FROM {table} WHERE {column} IN ({placeholders})"
    res = [row[0] for row in conn.execute(qry, ctx_set)]
    conn.close()
    return res


def querying_SQL_with_context(header: str, sql_col: str, sql_col_rev: str, ctx: list[str]):
    db = Path("/usr/src/app/pipeline/SQL_database") / f"{header}.sqlite"
    pos_5 = query_positions_sql(db, header, sql_col, ctx)
    pos_3 = query_positions_sql(db, header, sql_col_rev, ctx)
    return pos_5, pos_3


def get_mutations_table(header: str,
                        pos_5_sel: list[int],
                        pos_3_sel: list[int],
                        prob_table: pd.DataFrame) -> pd.DataFrame:
    db = Path("/usr/src/app/pipeline/SQL_database") / f"{header}.sqlite"
    conn = sqlite3.connect(db)

    def fetch_rows(pos_list):
        if not pos_list:
            return pd.DataFrame()
        placeholders = ",".join(map(str, pos_list))
        df = pd.read_sql_query(f"SELECT * FROM {header} WHERE ROWID IN ({placeholders})", conn)
        df["position"] = pos_list
        return df

    df5 = fetch_rows(pos_5_sel)
    df3 = fetch_rows(pos_3_sel)

    df5["Variant"] = df5["nucleotide"].apply(lambda n: generate_variant(n, prob_table))
    df3["Variant"] = df3["comp_nucleotide"].apply(lambda n: generate_variant(n, prob_table))
    df3["Variant"] = df3["Variant"].str.translate(str.maketrans("ACGT", "TGCA"))

    conn.close()
    return pd.concat([df5, df3], ignore_index=True)

###############################################################################
# V. MUTATED GENOME FRACTION OR MUTATION RATE?                               #
###############################################################################

def validate_mutation_mode(args: argparse.Namespace):
    if args.mut_genome_fraction is None and args.specific_mut_rate is None:
        sys.exit("One of --mut_genome_fraction or --specific_mut_rate must be provided")

    if args.mut_genome_fraction is not None:
        if not (0 < args.mut_genome_fraction <= 1):
            sys.exit("mut_genome_fraction must be >0 and ≤1")
        if args.fraction_positions is None:
            args.fraction_positions = 0.05
            print("Fraction of mutation not specified, using default of 0.05 (5%)")
        elif not (0 < args.fraction_positions <= 1):
            sys.exit("fraction_positions must be >0 and ≤1")
        print(f"Fraction of all genomes that will get a mutation is set to {args.mut_genome_fraction}")
        print(f"Fraction of positions that will be mutated is {args.fraction_positions}")

    if args.specific_mut_rate is not None:
        if not (0 <= args.specific_mut_rate <= 1):
            sys.exit("Invalid specific mutation rate")
        print(f"Specific mutation rate: {args.specific_mut_rate}")

###############################################################################
# MAIN PER‑SEQUENCE WORKER FUNCTION                                          #
###############################################################################

def worker(seq_rec: SeqRecord,
           copy_no: int,
           ctx: list[str] | None,
           ctx_len: int | None,
           sql_col: str,
           sql_col_rev: str,
           prob_table: pd.DataFrame,
           args: argparse.Namespace,
           worker_seed: int | None):

    if worker_seed is not None:
        random.seed(worker_seed)
        np.random.seed(worker_seed % (2**32 - 1))

    header = seq_rec.id
    fasta_len = len(seq_rec.seq)

    #######################  Mode A  #########################################
    if args.mut_genome_fraction is not None:
        # Determine candidate positions
        if ctx is not None:
            pos_5, pos_3 = querying_SQL_with_context(header, sql_col, sql_col_rev, ctx)
            candidate_positions = sorted(set(pos_5) | set(pos_3))
        else:
            candidate_positions = list(range(1, fasta_len + 1))
            pos_5 = pos_3 = []

        total_mut = round(len(candidate_positions) * args.fraction_positions)
        picked_positions = sorted(random.sample(candidate_positions, total_mut))

        # Split into 5' and 3' buckets if we have context information
        if ctx is not None:
            pos_5_sel = [p for p in picked_positions if p in pos_5]
            pos_3_sel = [p for p in picked_positions if p in pos_3]
        else:
            half = round(0.5 * len(picked_positions))
            pos_5_sel = random.sample(picked_positions, half)
            pos_3_sel = [p for p in picked_positions if p not in pos_5_sel]

        mutations_df = get_mutations_table(header, pos_5_sel, pos_3_sel, prob_table)
        mutations_df["pos_var"] = mutations_df["position"].astype(str) + mutations_df["Variant"]
        all_pos_vars = mutations_df["pos_var"].tolist()

        mg_number = round(copy_no * args.mut_genome_fraction)
        nmg_copies = copy_no - mg_number

        # Non‑mutated genomes (NMG)
        sample_rows = [
            {"genome_name": f"{header}_NMG", "variant_position": "No mutation", "copy_number": nmg_copies}
        ]

        # Mutated genomes (MG)
        mg_variants = {}
        for _ in range(mg_number):
            num_mut = random.choice(range(1, len(all_pos_vars) + 1))
            sel = sorted(random.sample(all_pos_vars, num_mut), key=lambda s: int("".join(filter(str.isdigit, s))))
            key = ",".join(sel)
            mg_variants[key] = mg_variants.get(key, 0) + 1

        for idx, (var_pos, count) in enumerate(mg_variants.items(), 1):
            sample_rows.append({
                "genome_name": f"{header}_MG{idx}",
                "variant_position": var_pos,
                "copy_number": count,
            })

        # Write outputs
        mutations_df.to_csv(f"{header}_mutations_table.csv", index=False)
        pd.DataFrame(sample_rows).to_csv(f"{header}_sample_table.csv", index=False)

    #######################  Mode B  #########################################
    else:
        lam = args.specific_mut_rate * fasta_len
        pois = np.random.poisson(lam, size=copy_no)
        has_mut = pois > 0
        n_mut_per_copy = pois[has_mut]

        if n_mut_per_copy.sum() == 0:
            print(f"No copies of genome {header} were mutated…")
            pd.DataFrame({
                "genome_name": [f"{header}_NMG"],
                "variant_position": ["No mutation"],
                "copy_number": [copy_no],
            }).to_csv(f"{header}_sample_table.csv", index=False)
            return

        # Prepare placeholder DataFrame for mutated genomes
        mg_rows = []
        if ctx is not None:
            pos_5, pos_3 = querying_SQL_with_context(header, sql_col, sql_col_rev, ctx)
            candidate_positions = sorted(set(pos_5) | set(pos_3))
        else:
            candidate_positions = list(range(1, fasta_len + 1))
            pos_5 = pos_3 = []

        all_selected_positions = set()
        for idx, n_mut in enumerate(n_mut_per_copy, 1):
            positions = sorted(random.sample(candidate_positions, n_mut))
            all_selected_positions.update(positions)
            mg_rows.append({
                "genome_name": f"{header}_MG{idx}",
                "nr_mutations": n_mut,
                "positions": ",".join(map(str, positions)),
            })

        # Build mutations table only once for the union of all positions
        if ctx is not None:
            pos_5_sel = [p for p in all_selected_positions if p in pos_5]
            pos_3_sel = [p for p in all_selected_positions if p in pos_3]
        else:
            pos_list = sorted(all_selected_positions)
            half = round(0.5 * len(pos_list))
            pos_5_sel = random.sample(pos_list, half)
            pos_3_sel = [p for p in pos_list if p not in pos_5_sel]

        mutations_df = get_mutations_table(header, pos_5_sel, pos_3_sel, prob_table)
        mutations_df["position"] = mutations_df["position"].astype(str)

        # Helper to map positions -> variant string
        pos2var = dict(zip(mutations_df["position"], mutations_df["Variant"]))
        for row in mg_rows:
            var_pos_list = [f"{p}{pos2var[str(p)]}" for p in map(int, row["positions"].split(","))]
            row["variant_position"] = ",".join(var_pos_list)
            row["copy_number"] = 1

        # Collapse identical variant strings
        df_mg = pd.DataFrame(mg_rows)
        df_sample_mg = (
            df_mg.groupby("variant_position", as_index=False)
            .agg(genome_name=("genome_name", "first"), copy_number=("copy_number", "sum"))
        )

        sample_df = pd.concat([
            pd.DataFrame({
                "genome_name": [f"{header}_NMG"],
                "variant_position": ["No mutation"],
                "copy_number": [copy_no - len(df_mg)],
            }),
            df_sample_mg,
        ], ignore_index=True)
        sample_df = sample_df[sample_df["copy_number"] > 0]

        # Write outputs
        mutations_df.to_csv(f"{header}_mutations_table.csv", index=False)
        sample_df.to_csv(f"{header}_sample_table.csv", index=False)

###############################################################################
# ENTRY POINT                                                                #
###############################################################################

def main():
    args = parse_args()

    # Change to output directory early
    args.output.mkdir(parents=True, exist_ok=True)
    os.chdir(args.output)
    print(f"Ok, we will use this {args.output} as the output folder, ok!")

    set_global_seed(args.set_seed)
    cores = configure_cores(args.cores)

    # Load input FASTA
    fasta_records = load_fasta(args.fasta_sequence)
    # Save header‑length CSV (fasta_lengths.csv)
    pd.DataFrame({
        "fasta_names": [rec.id for rec in fasta_records],
        "fasta_lengths": [len(rec.seq) for rec in fasta_records],
    }).to_csv("fasta_lengths.csv", index=False)
    print("fasta_length.csv - generated.")

    # Mutation context and probability table
    ctx = parse_mut_context(args.mutational_context)
    ctx_len = len(ctx[0]) if ctx else None
    sql_col, sql_col_rev = context_sql_columns(ctx_len)

    prob_table = build_variant_prob_table(args)

    # Copy number per FASTA header
    copies = determine_nr_copies(args, fasta_records)

    # Validate mutation mode parameters
    validate_mutation_mode(args)

    # Worker seeds if reproducible
    worker_seeds = None
    if args.set_seed is not None:
        worker_seeds = random.sample(range(1, 1_000_000_000), len(fasta_records))

    # Prepare multiprocessing pool
    with mp.Pool(processes=cores) as pool:
        pool.starmap(
            worker,
            zip(
                fasta_records,
                copies,
                repeat(ctx),
                repeat(ctx_len),
                repeat(sql_col),
                repeat(sql_col_rev),
                repeat(prob_table),
                repeat(args),
                worker_seeds or repeat(None),
            ),
        )

    print("Writing out the outputs!")
    print("Process finished")
    print("*_mutations_table.csv contains information about mutations inserted into mutated genomes")
    print("*_sample_table.csv contains the information about how many genomes got a mutation and which positions/nucleotide were mutated")
    print("Be aware… if the number of copies is too low as well as the mutation rate, it might happen that there are no mutated genomes at all, check and adjust the parameters!")


if __name__ == "__main__":
    main()
