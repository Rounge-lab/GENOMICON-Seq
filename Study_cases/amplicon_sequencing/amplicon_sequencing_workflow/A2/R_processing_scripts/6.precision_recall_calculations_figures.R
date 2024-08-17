library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(stringr)


key_table<-read.csv("A2_sample_key_index.csv", sep = "\t")

sample_tables<-read.csv("A2_sample_tables.csv")

#remove "_HPV16REF_sample_table" from csv_table columnin sample_tables
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

key_table_short<-key_table %>% select(sample_index,sample_table_id,pcr_table_id,sample_rep,tech_rep,number_copies,fragment_fraction)
key_table_short$sample_index<-gsub("-R|-F","", key_table_short$sample_index)
key_table_short<-key_table_short %>% group_by(sample_index) %>% distinct()

# add an identifier for the the original sample

key_table_short <- key_table_short %>%
  group_by(number_copies,fragment_fraction) %>%
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


### adding the reference nucleotide to the sample_mutation_count 

HPV16_table<-read.csv("HPV16_contexts.csv") %>% select(position, nucleotide,tri_context, tri_context_rev_comp)
sample_mutations_count$variant_position <-as.numeric(sample_mutations_count$variant_position)
sample_mutations_count<-sample_mutations_count %>% left_join(HPV16_table %>% select(position, nucleotide),by = c("variant_position"="position"))

APOBEC_targets<-HPV16_table %>% filter(tri_context == "TCA" |tri_context == "TCT"|tri_context == "TCG" | tri_context == "TCC"|
                                         tri_context_rev_comp =="TCA"| tri_context_rev_comp =="TCT"| tri_context_rev_comp =="TCG"|
                                         tri_context_rev_comp =="TCC") %>% select(position)

sample_mutations_count$mut_type<-ifelse((sample_mutations_count$nucleotide == "C" & sample_mutations_count$variant_nucleotide == "T") |
                                          (sample_mutations_count$nucleotide == "G" & sample_mutations_count$variant_nucleotide == "A"), "CtoT", "other")

sample_mutations_count$APOBEC_targets<-ifelse(sample_mutations_count$variant_position %in% APOBEC_targets$position, "AP_target", "other")
sample_mutations_count$AP_T_F<-ifelse(sample_mutations_count$mut_type == "CtoT" & sample_mutations_count$APOBEC_targets == "AP_target", "A3_Active","No")
## load the variants table 

all_variants<-read.csv("A2_all_variants.csv")

# set the cutoffs
read_cutoff<-seq(from = 5, to = 1000, by = 5)
read_cutoff<-c(1,read_cutoff)
#extract all the samples
all_samples<-unique(all_variants$sample_id2)

### make sample_Table to store the results

results_table_all_samples<-key_table_short %>% select(-sample_table_id,-pcr_table_id,-number_reads)

# Create a data frame for read_cutoff to facilitate crossing
cutoff_df <- tibble(read_cutoff = read_cutoff)

# Use crossing to expand df with each read_cutoff value
results_table_all_samples_all_cutoffs <- crossing(results_table_all_samples, cutoff_df)


results_table_all_samples_all_cutoffs$TP<-NA
results_table_all_samples_all_cutoffs$FP<-NA
results_table_all_samples_all_cutoffs$FN<-NA
results_table_all_samples_all_cutoffs$precision<-NA
results_table_all_samples_all_cutoffs$recall<-NA

#### loop will go for each sample and each cutoff


for (n in 1:length(all_samples)) {
  for (j in 1:length(read_cutoff)){
    current_cutoff<-read_cutoff[j]
    all_variants_filtered<-all_variants %>% filter(minor_base_count >= current_cutoff) %>% select(sample_id2,position,reference,minor_base)
    
    current_sample<-all_samples[n]
    current_sample_table_name<-key_table_short$sample_table_id[key_table_short$sample_index == current_sample]
    current_sample_table<-sample_mutations_count %>% filter(csv_table == current_sample_table_name) %>% select(-csv_table)
    
    current_variants_table<-all_variants_filtered %>% filter(sample_id2 == current_sample)  
    
    
    
    if (nrow(current_sample_table) == 0) {
      TP=0
      FP=nrow(current_variants_table)
      FN=0
      precision=0
      recall=0
      
    } else {
      
      # Perform a left join on the specified conditions
      joined_tables <- left_join(current_variants_table, current_sample_table %>% select(-occurrences), 
                                 by = c("position" = "variant_position", 
                                        "reference" = "nucleotide", 
                                        "minor_base" = "variant_nucleotide"))
      
      # Add the TP_FP column: 'TP' if the condition is met, 'FP' otherwise
      joined_tables <- joined_tables %>%
        mutate(TP_FP = ifelse(!is.na(mut_type), "TP", "FP"))
      
      # If you only want to retain the original columns of current_variants_table plus TP_FP
      result_table <- joined_tables %>%
        select(names(current_variants_table), TP_FP)
      
      TP=sum(result_table$TP_FP == "TP")
      FP=sum(result_table$TP_FP == "FP")
      FN=nrow(current_sample_table)-TP
      precision=TP/(TP+FP)
      recall=TP/(TP+FN)
    }
    
    results_table_all_samples_all_cutoffs$TP[results_table_all_samples_all_cutoffs$sample_index == current_sample & results_table_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-TP
    results_table_all_samples_all_cutoffs$FP[results_table_all_samples_all_cutoffs$sample_index == current_sample & results_table_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-FP
    results_table_all_samples_all_cutoffs$FN[results_table_all_samples_all_cutoffs$sample_index == current_sample & results_table_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-FN
    results_table_all_samples_all_cutoffs$precision[results_table_all_samples_all_cutoffs$sample_index == current_sample & results_table_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-precision
    results_table_all_samples_all_cutoffs$recall[results_table_all_samples_all_cutoffs$sample_index == current_sample & results_table_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-recall
    
  }
}

results_table_all_samples_all_cutoffs[is.na(results_table_all_samples_all_cutoffs)]<-0

## make summary per sample!

summary_all_mutations<-results_table_all_samples_all_cutoffs %>% group_by(original_sample_index, 
                                                                          number_copies,
                                                                          fragment_fraction,
                                                                          read_cutoff) %>% summarise(mean_precision = mean(precision),
                                                                                                 se_precision = sd(precision) / sqrt(n()),
                                                                                                 mean_recall = mean(recall),
                                                                                                 se_recall = sd(recall) / sqrt(n()),
                                                                                                 mean_TP = mean(TP),
                                                                                                 sd_TP = sd(TP),
                                                                                                 mean_FP = mean(FP),
                                                                                                 sd_FP = sd(FP),
                                                                                                 mean_FN = mean(FN),
                                                                                                 sd_FN = sd(FN),
                                                                                                 .groups = 'drop')


#### Plot 


summary_all_mutations$number_copies <- as.factor(summary_all_mutations$number_copies)
summary_all_mutations$fragment_fraction <- as.factor(summary_all_mutations$fragment_fraction)

# Calculate F1 Score and find the optimal cutoff
summary_all_mutations <- summary_all_mutations %>%
  group_by(original_sample_index) %>%
  mutate(F1 = 2 * (mean_precision * mean_recall) / (mean_precision + mean_recall),  
         optimal = F1 == max(F1, na.rm = TRUE)) %>%  
  ungroup()

summary_all_mutations<-summary_all_mutations %>%
  mutate(optimal = replace_na(optimal, FALSE))


summary_all_mutations<-summary_all_mutations %>%
  group_by(original_sample_index) %>%
  mutate(optimal = read_cutoff == min(read_cutoff[optimal])) %>%
  ungroup()


p <- ggplot(data = summary_all_mutations, aes(x = read_cutoff)) +
  geom_line(aes(y = mean_precision, color = "Precision"), size = 0.3) +
  geom_errorbar(aes(ymin = mean_precision - se_precision, ymax = mean_precision + se_precision, color = "Precision"), width = 0.01, size = 0.2) +
  geom_line(aes(y = mean_recall, color = "Recall"), size = 0.3) +
  geom_errorbar(aes(ymin = mean_recall - se_recall, ymax = mean_recall + se_recall, color = "Recall"), width = 0.01, size = 0.2) +
  geom_point(data = filter(summary_all_mutations, optimal), aes(y = mean_precision), color = "#57157e", size = 2) +
  geom_point(data = filter(summary_all_mutations, optimal), aes(y = mean_recall), color = "#ee5b5e", size = 2) +
  scale_y_continuous(name = "Precision",
                     sec.axis = sec_axis(~ ., name = "Recall")) +
  scale_color_manual(values = c("Precision" = "#57157e", "Recall" = "#ee5b5e")) +
  facet_wrap(~ original_sample_index, ncol = 4, nrow = 4, scales = "free_x") +
  labs(title = "Precision and Recall",
       x = "Read Cutoff",
       color = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none",
        strip.text = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())



ggsave("A2_precision_vs_recall_all_mut.jpg", p, height = 6, width = 6, units = "in")


### Narrowing the analysis to APOBEC3 signatures only


results_table_AP_all_samples_all_cutoffs<-key_table_short %>% select(-sample_table_id,-pcr_table_id,-number_reads)

# Create a data frame for read_cutoff to facilitate crossing
cutoff_df <- tibble(read_cutoff = read_cutoff)

# Use crossing to expand df with each read_cutoff value
results_table_AP_all_samples_all_cutoffs <- crossing(results_table_AP_all_samples_all_cutoffs, cutoff_df)


results_table_AP_all_samples_all_cutoffs$TP<-NA
results_table_AP_all_samples_all_cutoffs$FP<-NA
results_table_AP_all_samples_all_cutoffs$FN<-NA
results_table_AP_all_samples_all_cutoffs$precision<-NA
results_table_AP_all_samples_all_cutoffs$recall<-NA



APOBEC_variants<- all_variants %>%
  mutate(mut_type = case_when(
    (major_base == "C" & minor_base == "T") | (major_base == "G" & minor_base == "A") ~ "CtoT",
    TRUE ~ "other"
  ))

APOBEC_variants$APOBEC_targets<-ifelse(APOBEC_variants$position %in% APOBEC_targets$position, "AP_target", "other")
APOBEC_variants$AP_T_F<-ifelse(APOBEC_variants$APOBEC_targets == "AP_target" & APOBEC_variants$mut_type == "CtoT", "A3_active", "No")




for (n in 1:length(all_samples)) {
  for (j in 1:length(read_cutoff)){
    
    current_cutoff<-read_cutoff[j]
    all_variants_filtered<-APOBEC_variants %>% filter(minor_base_count >= current_cutoff) %>% select(sample_id2,position,reference,minor_base, AP_T_F)
    
    current_sample<-all_samples[n]
    current_sample_table_name<-key_table_short$sample_table_id[key_table_short$sample_index == current_sample]
    current_sample_table<-sample_mutations_count %>% 
      filter(csv_table == current_sample_table_name) %>% select(-csv_table) %>% filter(AP_T_F == "A3_Active")
    
    current_variants_table<-all_variants_filtered %>% filter(sample_id2 == current_sample) %>% filter(AP_T_F == "A3_active") %>% select(-AP_T_F)
    
    
    
    if (nrow(current_sample_table) == 0) {
      TP=0
      FP=nrow(current_variants_table)
      FN=0
      precision=0
      recall=0
      
    } else {
      
      # Perform a left join on the specified conditions
      joined_tables <- left_join(current_variants_table, current_sample_table %>% select(-occurrences), 
                                 by = c("position" = "variant_position", 
                                        "reference" = "nucleotide", 
                                        "minor_base" = "variant_nucleotide"))
      
      # Add the TP_FP column: 'TP' if the condition is met, 'FP' otherwise
      joined_tables <- joined_tables %>%
        mutate(TP_FP = ifelse(!is.na(AP_T_F), "TP", "FP"))
      
      # If you only want to retain the original columns of current_variants_table plus TP_FP
      result_table <- joined_tables %>%
        select(names(current_variants_table), TP_FP)
      
      TP=sum(result_table$TP_FP == "TP")
      FP=sum(result_table$TP_FP == "FP")
      FN=nrow(current_sample_table)-TP
      precision=TP/(TP+FP)
      recall=TP/(TP+FN)
    }
    
    results_table_AP_all_samples_all_cutoffs$TP[results_table_AP_all_samples_all_cutoffs$sample_index == current_sample & results_table_AP_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-TP
    results_table_AP_all_samples_all_cutoffs$FP[results_table_AP_all_samples_all_cutoffs$sample_index == current_sample & results_table_AP_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-FP
    results_table_AP_all_samples_all_cutoffs$FN[results_table_AP_all_samples_all_cutoffs$sample_index == current_sample & results_table_AP_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-FN
    results_table_AP_all_samples_all_cutoffs$precision[results_table_AP_all_samples_all_cutoffs$sample_index == current_sample & results_table_AP_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-precision
    results_table_AP_all_samples_all_cutoffs$recall[results_table_AP_all_samples_all_cutoffs$sample_index == current_sample & results_table_AP_all_samples_all_cutoffs$read_cutoff == current_cutoff]<-recall
    
  }
}

#########################
results_table_AP_all_samples_all_cutoffs[is.na(results_table_AP_all_samples_all_cutoffs)]<-0

## make summary per sample!

summary_APOBEC_mutations<-results_table_AP_all_samples_all_cutoffs %>% group_by(original_sample_index, 
                                                                          number_copies,
                                                                          fragment_fraction,
                                                                          read_cutoff) %>% summarise(mean_precision = mean(precision),
                                                                                                     se_precision = sd(precision) / sqrt(n()),
                                                                                                     mean_recall = mean(recall),
                                                                                                     se_recall = sd(recall) / sqrt(n()),
                                                                                                     mean_TP = mean(TP),
                                                                                                     sd_TP = sd(TP),
                                                                                                     mean_FP = mean(FP),
                                                                                                     sd_FP = sd(FP),
                                                                                                     mean_FN = mean(FN),
                                                                                                     sd_FN = sd(FN),
                                                                                                     .groups = 'drop')




#### Plot 

summary_APOBEC_mutations$number_copies <- as.factor(summary_APOBEC_mutations$number_copies)
summary_APOBEC_mutations$fragment_fraction <- as.factor(summary_APOBEC_mutations$fragment_fraction)

# filter out non-mutated samples

summary_APOBEC_mutations<-summary_APOBEC_mutations %>% filter(!original_sample_index %in% c(1,2,3,4))

# Calculate F1 Score and find the optimal cutoff
summary_APOBEC_mutations <- summary_APOBEC_mutations %>%
  group_by(original_sample_index) %>%
  mutate(F1 = 2 * (mean_precision * mean_recall) / (mean_precision + mean_recall),  
         optimal = F1 == max(F1, na.rm = TRUE)) %>%  
  ungroup()

summary_APOBEC_mutations<-summary_APOBEC_mutations %>%
  mutate(optimal = replace_na(optimal, FALSE))


summary_APOBEC_mutations<-summary_APOBEC_mutations %>%
  group_by(original_sample_index) %>%
  mutate(optimal = read_cutoff == min(read_cutoff[optimal])) %>%
  ungroup()

p <- ggplot(data = summary_APOBEC_mutations, aes(x = read_cutoff)) +
  geom_line(aes(y = mean_precision, color = "Precision"), size = 0.3) +
  geom_errorbar(aes(ymin = mean_precision - se_precision, ymax = mean_precision + se_precision, color = "Precision"), width = 0.01, size = 0.2) +
  geom_line(aes(y = mean_recall, color = "Recall"), size = 0.3) +
  geom_errorbar(aes(ymin = mean_recall - se_recall, ymax = mean_recall + se_recall, color = "Recall"), width = 0.01, size = 0.2) +
  geom_point(data = filter(summary_APOBEC_mutations, optimal), aes(y = mean_precision), color = "#57157e", size = 2) +
  geom_point(data = filter(summary_APOBEC_mutations, optimal), aes(y = mean_recall), color = "#ee5b5e", size = 2) +
  scale_y_continuous(name = "Precision",
                     sec.axis = sec_axis(~ ., name = "Recall")) +
  scale_color_manual(values = c("Precision" = "#57157e", "Recall" = "#ee5b5e")) +
  
  # Faceting
  facet_wrap(~ original_sample_index, ncol = 4, nrow = 4, scales = "free_x") +
  
  # Labels and themes
  labs(title = "Precision and Recall",
       x = "Read Cutoff",
       color = "Metric") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none",
        strip.text = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())



ggsave("A2_precision_vs_recall_APOBEC_mut.jpg", p, height = 6, width = 6, units = "in")



### save tables 

write.csv(x = results_table_all_samples_all_cutoffs, "A2_all_mutations_FP_TP_FN_per_tech_rep.csv", row.names = F)
write.csv(x = summary_all_mutations, "A2_all_mutations_FP_TP_FN_summary_per_sample.csv", row.names = F)
write.csv(x = results_table_AP_all_samples_all_cutoffs, "A2_APOBEC_FP_TP_FN_per_tech_rep.csv", row.names = F)
write.csv(x = summary_APOBEC_mutations, "A2_APOBEC_FP_TP_FN_summary_per_sample.csv", row.names = F)
