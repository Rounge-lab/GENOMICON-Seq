library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(stringr)
install.packages("devtools")
devtools::install_github("psyteachr/introdataviz")
library(introdataviz) # For geom_split_violin


key_table<-read.csv("A3_sample_key_index.csv", sep = "\t")

sample_tables<-read.csv("A3_sample_tables.csv")

#remove "_HPV16REF_sample_table" from csv_table column in sample_tables
sample_tables$csv_table<-gsub("_HPV16REF_sample_table", "", sample_tables$csv_table)
sample_tables$csv_table<-gsub("_config","",sample_tables$csv_table)

# add column for sample_tables identifier, fix the identifier to be identical to the table one!
key_table$sample_table_id<-sub("(-tech-replicate-.*)","", key_table$sample_id)
key_table$sample_table_id<-gsub("-","_", key_table$sample_table_id)
key_table$sample_table_id <- gsub("(\\d)e_(\\d)", "\\1e-\\2", key_table$sample_table_id, perl=TRUE)

# add column for pcr_table
key_table$pcr_table_id<-sub("(-HPV16-.*)","", key_table$sample_id)
key_table$pcr_table_id<-gsub("-","_", key_table$pcr_table_id)
key_table$pcr_table_id <- gsub("(\\d)e_(\\d)", "\\1e-\\2", key_table$pcr_table_id, perl=TRUE)


#make a shorter version key table

key_table_short<-key_table %>% select(sample_index,sample_table_id,pcr_table_id,sample_rep,tech_rep,copy_number,mutation_rate,polymerase_error_rate,number_reads)
key_table_short$sample_index<-gsub("-R|-F","", key_table_short$sample_index)
key_table_short<-key_table_short %>% group_by(sample_index) %>% distinct()

# add an identifier for the the original sample

key_table_short <- key_table_short %>%
  group_by(copy_number, mutation_rate, polymerase_error_rate, number_reads) %>%
  mutate(original_sample_index = cur_group_id()) %>%
  ungroup()

key_table_short$sample_rep<-NULL

# prepare sample_table

sample_mutations_count<- sample_tables %>% select(csv_table,Variant,position,occurrence)
sample_mutations_count$variant_position<-paste0(sample_mutations_count$position, sample_mutations_count$Variant)
sample_mutations_count<-sample_mutations_count %>% select(csv_table,variant_position, occurrence,Variant)  

colnames(sample_mutations_count)<-c("csv_table","variant_position","occurrences","variant_nucleotide")

sample_mutations_count <- sample_mutations_count %>%
  mutate(
    variant_nucleotide = sub(".*?(\\D)$", "\\1", variant_position),
    variant_position = sub("(\\d+).*", "\\1", variant_position)
  )

sample_mutations_count<-sample_mutations_count %>% select(csv_table, variant_position, variant_nucleotide, occurrences)

### Output how many unique mutations are there in total and what is the their mean occurence and SD - (occurence - number of genomes carrying a particular mutation)

summary_table <- sample_mutations_count %>%
  group_by(csv_table) %>%
  summarise(
    total_variant_positions = n_distinct(variant_position),
    mean_occurrences = mean(occurrences),
    sd_occurence = sd(occurrences)
  )

write.csv(summary_table, "A13_unique_mutations_overview.csv")


### Now load the variants - output from the script 4.alternative_allele_calling.R, for each sample, extract the true and the true positives
all_variants<-read.csv("A3_all_variants.csv")


## Try to filter the true positive 
all_variants$sample_id2 <- as.factor(all_variants$sample_id2)

all_variants <- all_variants %>%
  mutate(join_key = paste(position, minor_base))

mapped_ids <- key_table_short %>%
  select(sample_index, sample_table_id) %>%
  inner_join(all_variants, by = c("sample_index" = "sample_id2"))


sample_mutations_count <- sample_mutations_count %>%
  mutate(join_key = paste(variant_position, variant_nucleotide))

# Then perform the join on this composite key along with the sample_table_id mapping
true_variants <- mapped_ids %>%
  inner_join(sample_mutations_count, by = c("sample_table_id" = "csv_table", "join_key")) 

all_variants$join_key<-NULL

true_variants_adjusted <- true_variants %>%
  rename(sample_id2 = sample_index) %>%
  select(sample_id2, position, minor_base, minor_base_freq) # assuming minor_base_freq is not a factor for uniqueness

# Now, filter out false variants
false_variants <- anti_join(all_variants, true_variants_adjusted, 
                            by = c("sample_id2", "position", "minor_base", "minor_base_freq"))
false_variants<-false_variants %>% rename(sample_index = sample_id2)

# add the sample and replicate identifier to the true and false variants

true_variants <- true_variants %>%
  left_join(key_table_short %>% select(sample_index, original_sample_index, tech_rep), by = c("sample_index"))


false_variants<-false_variants %>% 
  left_join(key_table_short %>% select(sample_index, original_sample_index, tech_rep), by = c("sample_index"))

#### extract sample ids with parameters

sample_indexes_table<-key_table_short %>% 
  select(original_sample_index,copy_number,mutation_rate,polymerase_error_rate,number_reads) %>% 
  distinct() %>% arrange(original_sample_index)

### load context table containing the info about the trinucleotide context for each positions
HPV16_context<-read.csv("HPV16_contexts.csv")
HPV16_context <-HPV16_context %>% select(position,nucleotide,comp_nucleotide,tri_context,tri_context_rev_comp)


## adding additional column true_false

true_variants<-true_variants%>% mutate(true_false = "true")
false_variants<-false_variants%>% mutate(true_false = "false")

true_variants<- true_variants %>% select(reference,sample_index,chr,position,A,G,C,T,CumCoverage,
                                         major_base,major_base_count,minor_base,minor_base_count,
                                         coverage,minor_base_freq,original_sample_index,tech_rep,
                                         true_false)
## combine true and false variants

true_and_false<-rbind(false_variants,true_variants)

true_and_false<-true_and_false %>% select(original_sample_index,tech_rep,position, minor_base_count, true_false)

false_color <- "#fa7f5e" 
true_color <- "#1fa287"  

gg <- ggplot(true_and_false, aes(x = original_sample_index, y = minor_base_count, fill = true_false)) +
  introdataviz::geom_split_violin(trim = TRUE, alpha = 0.4, scale = "area") +
  geom_boxplot(width = 0.13, position = position_dodge(width = 0.2), alpha = 1, fatten = NULL, outlier.shape = NA) +
  scale_y_log10(limits = c(1, 10000)) +  
  scale_fill_manual(values = c("false" = false_color, "true" = true_color)) +
  labs(y = "Minor Base Count", x = "Original Sample Index") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = unit(c(0,0,0,0), "lines")
  )


ggsave("A3_mutations_vs_noise_split_violin_all_mutation.jpg", gg, width = 10, height = 3, units = "in")

summary_table <- true_and_false %>%
  group_by(original_sample_index, true_false) %>%
  summarise(
    count_observations = n()/5,  # mean count of observations for 5 replacates
    mean_read_count = mean(minor_base_count, na.rm = TRUE),
    min_read_count = min(minor_base_count, na.rm = TRUE),
    max_read_count = max(minor_base_count, na.rm = TRUE),
    sd_read_count = sd(minor_base_count, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  # Pivot to get true/false metrics into separate columns
  pivot_wider(names_from = true_false, values_from = c(count_observations, mean_read_count, min_read_count, max_read_count, sd_read_count)) %>%
  # Rename for clarity
  rename_with(.cols = everything(), .fn = ~ str_replace(., "_observations", ""))



write.csv(summary_table,"A3_summary_read_counts_counts_all_mutations.csv",)


