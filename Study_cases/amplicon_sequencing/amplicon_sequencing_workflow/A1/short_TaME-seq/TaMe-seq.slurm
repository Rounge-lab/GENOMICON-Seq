#!/bin/bash
#SBATCH --account=ec34
#SBATCH --job-name=Tameseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50            
#SBATCH --mem-per-cpu=6G                    
#SBATCH --time=06:00:00

module purge

source ../TaME-seq_env/bin/activate

snakemake -j 48 -p -s ./Snakefile_tameseq_short_mod_cutadapt --config directory=run_7
