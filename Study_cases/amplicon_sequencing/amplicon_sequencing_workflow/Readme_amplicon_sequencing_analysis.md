### Running the ampliseq A1-A3 study case

Each study case has its own folder, and although the scripts are very similar, they are not 100% identical as study cases each had their own variable parameters used to generate simulated data.
Here, we present how we run the GENOMICON-Seq on EduCloud HPC, and how we analyzed the generated samples in this environment, running scripts on a personal computer or other HPV would require adjustments.

1. script `1.making_configs.sh` 
After setting up the GENOMICON-Seq for amplicon sequencing, place the `1.making_configs.sh` into the folder, and run `./1.making_configs.sh`
The script requires an array of values for each customizable parameter to be set up. It uses the one a config_ampliseq.yml as the template (the parameters in the template would then be used in all the versions of the made configs). 
The script then creates a set of config files with the specified parameters - one config corresponds to one sample that will be generated with the given parameters.

2. script `2.is_ampliseq_many_configs.slurm` and `3.sumitting_jobs.sh` 
script `2.is_ampliseq_many_configs.slurm` is the template for the `3.sumitting_jobs.sh` script that submits the job to HPC, in this case the scripts were made to be compatible with the EduCloud HPC we used `https://www.uio.no/english/services/it/research/platforms/edu-research/`. `3.sumitting_jobs.sh` combines the `2.is_ampliseq_many_configs.slurm` script and each config file, creates a separate folder for each GENOMICON-Seq run with the symbolic link to all the input files, and submits the job.

3. After the jobs are done, remove the created symbolic links in each working folder.

4. Copy and past config files to the corresponding working folder

5. script `4.renaming_and_copying_FASTQ_files.sh` 
The script searches for all the produced FASTQ files, in all the folders, uses the config file name to rename the FASTQ file so that they will contain the parameters and their values that were varying in this study. The file names has particular structure to be compatible with the TaME-seq pipeline

6. script `5.extracting_mutations_overview_table.sh`
The script extracts the csv table (_inserted_mutations_overview.csv) from every working GENOMICON-seq folder. The table contains the overview of the inserted mutations in each generated sample! Make sure to extract the tables in their own folder.

7. script `6.combining_mutations_overview_table.sh`
The script can be used to combine all the extracted tables (inserted_mutations_overview.csv) into one with the addition of the original table name as column so that each mutation can be traced to the correct sample they originate from

8. script `7.extracting_mutations_counts_table.sh`
The script extracts the csv table (_mutations_counts.csv) from every working GENOMICON-seq folder. The table contains the overview of the number of mutations each genome copy got during sample generation! Make sure to extract the tables in their own folder.

9. script `8.combining_mutations_counts_table.sh`
The script can be used to combine all the extracted tables (mutations_counts.csv) into one with the addition of the original table name as a new column.

10. script `9.extracting_pcr_mutations_table.sh` and `10.merging_pcr_mutations_table.sh` search for the `PCR_*-polymerase_error_mutations.csv` tables (depending on the number of PCR reaction there will be one table with the info on which mutation were mutated during PCR).

11. script `11.combining_merged_pcr.sh` 
The script will combine the info from two PCR reactions into one! - MAYBE REMOVE

12. Transfer the extracted and renamed FASTQ files into the TaME-seq input data folder and run the TaME-seq (files and folders of shorter version of the TaME-seq as well as the slurm script for job submission that is not performing the analysis of the HPV-integration is in the `short_TaME-seq` folder)

13. After the TaME-seq run is finished, the output files (*.fastq_stats.csv, *.chromosome.csv, *.coverage.csv, *.stats.csv) are further processed. 

14. script `12.create_key_table.csv`
The script is used with the output table *.fastq_stats.csv and outputs the `sample_key_index.csv` (specifiable as the script argument). Table will contain the original sample names, short sample names (), also for each produced sample (in this case each sample corresponds to the sequenced technical replicate of a generated in-silico sample) a parameter overview (variable parameter relevant for the study case) will be made.

15. script `renaming_sample_id.sh`
The script is is applied on all the TaME-seq outputs with the `sample_key_index.csv` to change the sample name into the shorter format stored in the `sample_key_index.csv`

16. In the next step we R-scripts in the `R_processing_scripts` folder will be used for further analysis.The scripts have to be run from Rstudio!

17. script `edit_overhangs.R`
TaME-seq pipeline maps the reads to the reference HPV which is elongated so that the sequence ends are extended by 1000 bp coming from the opposite ends. Since HPV is circular this is necessary to ensure that the reads would map correctly. The script actually combines the mapped nucleotides on the extended ends so that the correct ends! 

18. script `2.compiled_stats.R`
Uses output from the TaME-seq run and combine them into a csv table with the information on the mapping statistics of each sample replicate (for example, number of reads, number of mapped reads, number of unmapped reads, etc.)

19. script `3.counts_coverage.R`
Uses output from the TaME-seq run and extracts the sequencing depth statistics, such as mean/median sequencing depth across the HPV genome, and the % of the genome covered.

20. script `4.alternative_allele_calling.R`
The script preforms the filtering step of the coverage.csv table with the information about how many reads support the most common allele, and alternative allele in all positions on the reference genome. The filtering step is the removal of all the alternative alleles which mean reads quality is bellow 30. Filtered alternative alleles will be used in the final script to perform iSNV calling at different read-count cutoffs (number of reads supporting the alternative).

21. script `5.read_count_plot_true_vs_noise.R`
Script performs the analysis of the read counts supporting the mutations which signature matches the true mutation and mutations not matching these signatures categorized as noise. The main output is the figure showing the difference in the read count distributions supporting noise or generated mutations.

22. script `6.precision_recall_calculations_figures.R`
Script classifies the mutations into true positives, true negatives and false positives, and calculates the precision and recall for all the read-count cutoffs and identifies the optimal count cutoff, and outputs the plots and summary tables.

23. script `7.prec_rec_plot_optimal_read_count.R`
Script creates a plot based on the summary table made in the previous step by filtering the maximized precision and recall at optimal read-count cutoffs




