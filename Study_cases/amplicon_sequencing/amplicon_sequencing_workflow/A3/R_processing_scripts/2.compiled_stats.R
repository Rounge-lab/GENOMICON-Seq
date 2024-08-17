  
#think that reshape two needs to be installed :)
library("reshape2")
library("tidyr")
library("dplyr")
library(cowplot)

### Upload files
# .stat.csv table
table <-read.csv("Raw/A3.HPV16.stats.csv", header = TRUE, sep = "\t")

#chromosome.csv table
table_chr <- read.csv("Raw/A3.HPV16.chromosome.csv", header = TRUE, sep = "\t")

#fastq_stat.csv
fastq <- read.csv("Raw/A3.fastq_stats.csv",sep = "\t")


### Create output file

stats <- spread(table, stats, value)

stats$reference <- stats$strain


table_chr$target <- table_chr$strain

table_chr$target <- ifelse(table_chr$strain == table_chr$target,"target", "no")
table_chr$target[is.na(table_chr$target)]  <- "no"
table_chr$mapping <- ifelse(grepl("chr", table_chr$chr), "human", "hpv")

stats_chr <- table_chr %>%
  group_by(sample_id, strain, strand, mapper) %>%
  summarise(
    reads_mapped_human = sum(count[mapping=="human"]),
    reads_mapped_hpv = sum(count[mapping=="hpv"]),
    reads_mapped_target = sum(count[mapping=="hpv"]),
    reads_mapped_target = sum(count[target=="target"])
  )




fastq$raw_reads <- fastq$raw_reads*2
fastq$trimmed_reads <- fastq$trimmed_reads*2
fastq$X<-NULL

#takes the relevant variables from stats table and makes a stats_filt table
stats_filt <- stats[c("sample_id","reference","mapper","reads mapped","reads unmapped","insert size average","insert size standard deviation")]
#fix the names of the columns
colnames(stats_filt)[colnames(stats_filt) == 'reads mapped'] <- 'reads_mapped'
colnames(stats_filt)[colnames(stats_filt) == 'reads unmapped'] <- 'reads_unmapped'
colnames(stats_filt)[colnames(stats_filt) == 'insert size average'] <- 'insert_size'
colnames(stats_filt)[colnames(stats_filt) == 'insert size standard deviation'] <- 'insert_size_SD'

stats_merged <- merge(merge(fastq, stats_filt, by=c("sample_id"), all=T), stats_chr, by=c("sample_id","mapper"), all=T)

stats_merged <- stats_merged[,c("sample_id","strain","strand","reference","mapper","raw_reads","trimmed_reads",
                                "raw_length","r1_trimmed_length","r2_trimmed_length","insert_size","insert_size_SD",
                                "reads_mapped","reads_unmapped")]

write.csv(file="A3_compiled_stats.csv", stats_merged, row.names = FALSE)


