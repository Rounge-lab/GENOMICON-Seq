args <- commandArgs(trailingOnly = TRUE)

# Check if the number of arguments is as expected
if(length(args) != 4) {
  stop("Usage: Rscript script_name.R <coverage_file> <exon_file> <key_table_file> <output_file>", call. = FALSE)
}

# Assign arguments to variables
coverage_file <- args[1]
exon_file <- args[2]
key_table_file <- args[3]
output_file <- args[4]

library("ggplot2")
library("reshape2")
library("splitstackshape")
library("tidyr")
library("dplyr")
library(purrr)
library(data.table)



coverage <- read.csv(coverage_file, header = T, sep = "\t")  %>% select(chr,position,reference,coverage,sample_id)
# read huamn exom files

exons<-read.csv(exon_file, sep = "\t", header = T, col.names = c("chr", "start", "end"))
setDT(exons)
exons[, exon_length := (end - start + 1)]
exons[, exon_id := .I]

setDT(coverage)
coverage_sum <- coverage[exons, 
                         on = .(chr, position >= start, position <= end),
                         .(sample_id, exon_id = i.exon_id, coverage),
                         by = .EACHI]

# Aggregate to sum up coverage for each exon_id and sample_id combination
coverage_sum_aggregated <- coverage_sum[, .(total_coverage = sum(coverage)), by = .(sample_id, exon_id)]

# Merge aggregated coverage data with exons to get the exon_length
coverage_with_length <- merge(coverage_sum_aggregated, exons[, .(exon_id, exon_length)], by = "exon_id")

# Calculate mean coverage by dividing total coverage by exon_length for each sample for each exon
coverage_with_length[, mean_coverage := total_coverage / exon_length]

# Select the necessary columns
final_result <- coverage_with_length[, .(sample_id, exon_id, mean_coverage, exon_length)]
final_result$mean_coverage[is.na(final_result$mean_coverage)] <- 0



final_result <- as.data.table(final_result)

# Identify all unique sample_ids and exon_ids
unique_sample_ids <- unique(na.omit(final_result$sample_id))
unique_exon_ids <- unique(final_result$exon_id)

template_dt <- CJ(sample_id = unique_sample_ids, exon_id = unique_exon_ids)
template_dt<-merge(template_dt, exons %>% select(exon_id,exon_length), by=c("exon_id"))

expanded_dt <- merge(template_dt, final_result, by = c("sample_id", "exon_id", "exon_length"), all.x = TRUE)
expanded_dt[is.na(mean_coverage), mean_coverage := 0]

expanded_dt<-expanded_dt %>%
  mutate(length_category_when_not_covered = case_when(
    mean_coverage == 0 & exon_length > 1 & exon_length <= 50 ~ "1to50",
    mean_coverage == 0 & exon_length > 50 & exon_length <= 100 ~ "51to100",
    mean_coverage == 0 & exon_length > 100 & exon_length <= 200 ~ "101to200",
    mean_coverage == 0 & exon_length > 200 ~ ">200",
    TRUE ~ "covered_or_length_not_in_criteria"
  ))

total_number_exons<-nrow(exons)

final_result_summary <- expanded_dt %>%
  group_by(sample_id) %>%
  summarise(
    count_greater_than_10 = (sum(mean_coverage >= 10, na.rm = TRUE)/total_number_exons)*100,
    count_greater_than_50 = (sum(mean_coverage >= 50, na.rm = TRUE)/total_number_exons)*100,
    count_greater_than_100 = (sum(mean_coverage >= 100, na.rm = TRUE)/total_number_exons)*100,
    count_greater_than_200 = (sum(mean_coverage >= 200, na.rm = TRUE)/total_number_exons)*100,
    count_greater_than_250 = (sum(mean_coverage >= 250, na.rm = TRUE)/total_number_exons)*100,
    count_greater_than_300 = (sum(mean_coverage >= 300, na.rm = TRUE)/total_number_exons)*100,
    not_covered_count = (sum(mean_coverage == 0, na.rm = TRUE)/total_number_exons)*100,
    not_covered_length_1to50 = (sum(length_category_when_not_covered == "1to50", na.rm = TRUE)/total_number_exons)*100,
    not_covered_length_51to100 = (sum(length_category_when_not_covered == "51to100", na.rm = TRUE)/total_number_exons)*100,
    not_covered_length_101to200 = (sum(length_category_when_not_covered == "101to200", na.rm = TRUE)/total_number_exons)*100,
    not_covered_length_bt_200 = (sum(length_category_when_not_covered == ">200", na.rm = TRUE)/total_number_exons)*100
  )

## add a key table

key_table<-read.csv(key_table_file, sep="\t")

final_result_summary_merged_key<-left_join(final_result_summary, key_table %>% select(-sample_id,-sample_number), by=c("sample_id" = "sample_index"))

write.csv(final_result_summary_merged_key,output_file, row.names = F)
