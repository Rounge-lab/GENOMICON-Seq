#!/bin/bash
#SBATCH --account=ec34
#SBATCH --job-name=wes-Tameseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30            
#SBATCH --mem-per-cpu=16G                    
#SBATCH --time=48:00:00

module purge

source ../wes_tame_seq/bin/activate

snakemake -j 28 -p -s ./Snakefile_TaMEseq_wes --config directory=Run_17
