### Running and analysis of WES study cases

This the list of scripts used to process all files in W1-4 WES study cases. 
 
1. script `1.making_configs.sh` 
After setting up the GENOMICON-Seq for WES, place the `1.making_configs.sh` into the folder, and run `./1.making_configs.sh`
The script requires setting up values for each parameter. It uses the one a config_wes.yml as the template (the parameters in the template would then be used in all the versions of the made configs). 
The script then creates a set of config files with the specified parameters - one config corresponds to one sample that will be generated with the given parameters. Remember to make one mutation-free (no generated mutations) sample with the same parameters (number of copies, number of reads, etc but no mutation mode) because it will be needed during the Mutect2 run. This sample will be a control (non-tumor) sample. 

2. script `2.wes_many_configs.slurm` and `3.sumitting_jobs.sh` 
script `2.wes_many_configs.slurm` is the template for the `3.sumitting_jobs.sh` script that submits the job to HPC, in this case the scripts were made to be compatible with the EduCloud HPC we used `https://www.uio.no/english/services/it/research/platforms/edu-research/`. `3.sumitting_jobs.sh` combines the `2.wes_many_configs.slurm` script and each config file, creates a separate folder for each GENOMICON-Seq run with the symbolic link to all the input files, and summits the job.

3. After the jobs are done, remove the created symbolic links in each working folder.

4. Copy and past config files to the corresponding working folder

5. script `4.renaming_and_copying_FASTQ_files.sh` 
The script searches for all the produced FASTQ files, in all the folders, uses the config file name to rename the FASTQ file so that they will contain the parameters and their values that varied. The file names has particular structure to be compatible with the TaME-seq pipeline. Although TaME-seq pipeline was originally made for HPV processing, it has been modified to perform the mapping and extract the information of the mapping statistics. In addition, Bowtie2 mapper replaced the HISAT2.

6. script `5.extracting_mutations_overview_table.sh`
The script extracts the csv table (_inserted_mutations_overview.csv) from every working GENOMICON-seq folder. The table contains the overview of the inserted mutations in each generated sample! Make sure to extract the tables in their own folder.

7. script `6.combining_mutations_overview_table.sh`
The script can be used to combine all the extracted tables (inserted_mutations_overview.csv) into one with the addition of the original table name as column so that each mutation can be traced to the correct sample they originate from

8. script `7.extracting_mutations_counts_table.sh`
The script extracts the csv table (_mutations_counts.csv) from every working GENOMICON-seq folder. The table contains the overview of the number of mutations each genome copy got during sample generation! Make sure to extract the tables in their own folder.

9. script `8.combining_mutations_counts_table.sh`
The script can be used to combine all the extracted tables (mutations_counts.csv) into one with the addition of the original table name as a new column.

12. Transfer the extracted and renamed FASTQ files into the TaME-seq input data folder and run the modified TaME-seq pipeline (files for modified version of the TaME-seq as well as the slurm script for job submission is in the `mapping` folder)

13. After the mapping run is finished, the output files (*.fastq_stats.csv, *.chromosome.csv, *.coverage.csv, *.stats.csv) are further processed.

14. script `create_sample_key_table_wes.sh` (`mapping` folder)
Use script with the TaME-seq-generated table `*.fastq_stats.csv` to generate the short sample indexes and place them in the `key_table.csv`. Since the original sample names contain the information about the varying parameters, these will be extracted and placed in their own columns. 

15. scripts `renaming_sample_id.sh` and `renaming_coverage.py` (`mapping` folder)
`renaming_sample_id.sh` is used to replace the original sample names with the sample names in the in the TaME-Seq outputs *.fastq_stats.csv, *.chromosome.csv, *.stats.csv. Since *.coverage.csv table is much larger, `renaming_coverage.py` is used to to rename the sample names in this table (python makes it much faster).

16. R-scripts for quick extraction of the mapping statistics and coverage (`mapping/R_processing` folder)
`1compiled_stats_wes.R` needs to be run from Rstudio, it will output the basic mapping statistics.

`2counts_coverage_clean_wes.R` can be run outside of Rstudio, requires the BED file with info on exonic regions targeted in the analysis. Output the table with the coverage status over the whole exome (number of covered exons and the sequencing depth of the exons). `wes-R_Rscript_2.slurm` can be used to submit the `2counts_coverage_clean_wes.R` to the HPC so that the processing time of `2counts_coverage_clean_wes.R` would be faster.

17. Mutect2 run (`Mutect2_and_res_analysis` folder)
Script used for Mutect2 run are in the `Mutect2_run_and_analysis` folder. 
Download the Mutect2 required files:

gnomad
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

PON
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi

vcf.gz files can be filtered to include only hg38 chromosomes of interest, but then they need to be indexed (.tbi) with `tabix -p vcf file.vcf.gz`. In our case, the filtered chr1 files and their index files was placed in the folder necessary_files.

BAM-files (TaME-seq pipeline using Bowtie2) were placed in `bams` folder. Edit the key_table.csv so that it will contain all the necessary info for running the bash script that runs the Mutect2 container. Example is in the folder `Mutect2_run_and_analysis`.

Table contains the original BAM file name without the suffix, the suffix added to all file names when bam is created, tagged bam suffix (AddOrReplaceReadGroups continue reading), sample index corresponding to renamed samples, then varying parameters and their values, original sample number and its corresponding replicates, and in the end the the name that the sample was tagged with during AddOrReplaceReadGroups (this part actually distinguish the tumor and control samples and instructs the script running the Mutect2 which tumor sample should be compared with which normal sample).

18. script `mutect2_run.sh` and `mutect2.slurm`
`mutect2_run.sh` script first uses AddOrReplaceReadGroups to tag the files using the tag from the last column of the key_table.csv, then identifies pairs of samples (control-tumor) and run the Mutect2. The script is parallelized so that it can run Mutect2 call on different pairs of control-tumor simultaneously. 

`mutect2.slurm` script sends the job to the HPC.

19. script `filter_vcf_convert_to_csv.sh`
The script uses FilterMutectCalls (basic filtering) and converts the filtered vcf file into a csv table.

20. Script `vcf_csv_processing_and_figures.R` (`Mutect2_run_and_analysis/R_analysis` folder) - used for W1-3 study cases
Final processing R script (to be run from Rstudio), requires the merged mutations_overview.csv table (listing all the mutations in the original generated samples). It identifes the TP, FP, and FN for each sequenced replicate, and calculates the precision and recall. The script also creates the plots with the precision and recall changes over different parameters used to generate samples in GENOMICON-seq.   

21. Script `vcf_csv_processing_and_figures_w4.R`
For the study case W4, this script identifies the TPs, and plots their tri-mutational context. It also plots the occurrence of the tri-mutational context of the mutations in the generated sample, and the original SBS mutational signature based on the info from the table from the COSMIC database.  