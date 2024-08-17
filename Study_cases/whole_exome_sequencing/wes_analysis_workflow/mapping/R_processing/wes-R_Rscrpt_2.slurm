#!/bin/bash
#SBATCH --account=ec34
#SBATCH --job-name=wes-R
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2            
#SBATCH --mem-per-cpu=48G                    
#SBATCH --time=2:30:00

module purge
module load R/4.2.1-foss-2022a

Rscript ./scripts/3counts_coverage_clean_v3_wes_outside_args.R ./Raw/Run_11.chr1.coverage.csv ./chr1.bed ./key_table.csv ./R_res/coverage_stats_per_sample_per_rep_RUN_11.csv

