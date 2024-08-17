  
#think that reshape two needs to be installed :)
library("reshape2")
library("tidyr")
library("dplyr")
library(cowplot)

### Upload files
#FIX THE AUTOMATIC INPUT!
# .stat.csv table
table <-read.csv("Raw/Run_11.chr1.stats.csv", header = TRUE, sep = "\t")

#chromosome.csv table
table_chr <- read.csv("Raw/Run_11.chr1.chromosome.csv", header = TRUE, sep = "\t")

#fastq_stat.csv
fastq <- read.csv("Raw/Run_11.fastq_stats.csv",sep = "\t")


### Create output file

stats <- spread(table, stats, value)

stats$strain<-NULL

fastq$raw_reads <- fastq$raw_reads*2
fastq$trimmed_reads <- fastq$trimmed_reads*2
fastq$X<-NULL

#takes the relevant variables from stats table and makes a stats_filt table
stats_filt <- stats[c("sample_id","mapper","reads mapped","reads unmapped","insert size average","insert size standard deviation")]
#fix the names of the columns
colnames(stats_filt)[colnames(stats_filt) == 'reads mapped'] <- 'reads_mapped'
colnames(stats_filt)[colnames(stats_filt) == 'reads unmapped'] <- 'reads_unmapped'
colnames(stats_filt)[colnames(stats_filt) == 'insert size average'] <- 'insert_size'
colnames(stats_filt)[colnames(stats_filt) == 'insert size standard deviation'] <- 'insert_size_SD'

stats_merged <- merge(fastq, stats_filt, by=c("sample_id"))

stats_merged <- stats_merged[,c("sample_id","raw_reads","trimmed_reads",
                                "raw_length","r1_trimmed_length","r2_trimmed_length","insert_size","insert_size_SD",
                                "reads_mapped","reads_unmapped")]
#combine with key table!
key_table<-read.csv("key_table.csv", sep = "\t")

stats_merged_key<-left_join(stats_merged, key_table %>% dplyr::select(-sample_id), by=c("sample_id" = "sample_index"))
stats_merged_key<-stats_merged_key %>% dplyr::select(sample_id,original_sample_index,tech_rep,copy_number, mutation_rate,read_number,raw_reads,reads_mapped,reads_unmapped)


### now the same table as a summary per sample (mean, sd) for raw read, mapped reads and unmapped reads

stats_merged_key_per_original_sample<-stats_merged_key%>%
  group_by(original_sample_index, read_number) %>%
  summarise(
    copy_number = first(copy_number),
    mutation_rate = first(mutation_rate),
    read_number = first(read_number),
    mean_raw_reads = mean(raw_reads),
    sd_raw_reads = sd(raw_reads),
    mean_reads_mapped = mean(reads_mapped),
    sd_reads_mapped = sd(reads_mapped),
    mean_reads_unmapped = mean(reads_unmapped),
    sd_reads_unmapped = sd(reads_unmapped)
    
  )

stats_merged_key_per_original_sample_arranged<-stats_merged_key_per_original_sample%>%
  arrange(read_number, copy_number, mutation_rate)
stats_merged_key_per_original_sample_arranged$original_sample_index<-1:36


stats_merged_key_per_original_sample_arranged<-stats_merged_key_per_original_sample_arranged %>% select(original_sample_index,copy_number,mutation_rate,mean_raw_reads, sd_raw_reads, mean_reads_unmapped, sd_reads_unmapped)


write.csv(stats_merged_key_per_original_sample_arranged, "R_res/Run11_comp_stats_summary_per_table.csv", row.names = F)


coverage_stats<-read.csv("R_res/coverage_stats_per_sample_per_rep_RUN_11.csv")

coverage_stats<-coverage_stats %>% select(original_sample_index,tech_rep,read_number,copy_number,mutation_rate,count_greater_than_10,
                                          count_greater_than_50,count_greater_than_100,count_greater_than_200,count_greater_than_250,count_greater_than_300,
                                          not_covered_count, not_covered_length_1to50,not_covered_length_51to100,not_covered_length_101to200, not_covered_length_bt_200)

coverage_stats_summary<-coverage_stats %>% group_by(original_sample_index, read_number, copy_number,mutation_rate) %>% 
  summarise(
    mean_count_greater10 = mean(count_greater_than_10),
    sd_count_greater_10 = sd(count_greater_than_10),
    mean_count_greater50 = mean(count_greater_than_50),
    sd_count_greater_50 = sd(count_greater_than_50),
    mean_count_greater100 = mean(count_greater_than_100),
    sd_count_greater_100 = sd(count_greater_than_100),
    mean_count_greater200 = mean(count_greater_than_200),
    sd_count_greater_200 = sd(count_greater_than_200),
    mean_count_greater300 = mean(count_greater_than_300),
    sd_count_greater_300 = sd(count_greater_than_300),
    mean_not_covered = mean(not_covered_count),
    sd_not_covered = sd(not_covered_count),
    mean_not_covered_1to50 = mean(not_covered_length_1to50),
    sd_not_covered_1to50 = sd(not_covered_length_1to50),
    mean_not_covered_51to100 = mean(not_covered_length_51to100),
    sd_not_covered_51to100 = sd(not_covered_length_51to100),
    mean_not_covered_101to200 = mean(not_covered_length_101to200),
    sd_not_covered_101to200 = sd(not_covered_length_101to200),
    mean_not_covered_bt200 = mean(not_covered_length_bt_200),
    sd_not_covered_bt_200 = sd(not_covered_length_bt_200)
  )

coverage_stats_summary<-coverage_stats_summary %>%
  arrange(read_number, copy_number, mutation_rate)
coverage_stats_summary$original_sample_index<-1:36

write.csv(coverage_stats_summary, "R_res/Run11_coverage_stats_per_sample.csv", row.names = F)
