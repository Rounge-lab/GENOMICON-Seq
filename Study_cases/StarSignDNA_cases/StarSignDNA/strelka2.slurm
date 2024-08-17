#!/bin/bash
#SBATCH --account=ec34
#SBATCH --job-name=strelka2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22            
#SBATCH --mem-per-cpu=16G                    
#SBATCH --time=2:00:00

module purge
module load Miniconda3/22.11.1-1
conda activate /fp/homes01/u01/ec-milas/.conda/envs/strelka2_env

./strelka2.sh ./normals/sample_1_mut_r_0_perfect_r_35000000_c_1000_sorted.bam ./perfect_bams ./ref_genome/human_genome.fasta ./ref_genome/xgen-exome-research-panel-v2-probes-hg38_short_all.bed.gz 20

