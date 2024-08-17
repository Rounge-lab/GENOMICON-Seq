#!/bin/bash
#SBATCH --account=ec34
#SBATCH --job-name=MUTECT2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20            
#SBATCH --mem-per-cpu=24G                    
#SBATCH --time=64:00:00

module purge
module load Miniconda3/22.11.1-1
conda activate /fp/homes01/u01/ec-milas/.conda/envs/mutect2_add_stuff # conda env containing smatools necessary for indexing of the newly tagged bams

./mutect2_run.sh bams output output/bams_tagged output/vcf_output 18 key_table_mod_order_35M_3500.csv

