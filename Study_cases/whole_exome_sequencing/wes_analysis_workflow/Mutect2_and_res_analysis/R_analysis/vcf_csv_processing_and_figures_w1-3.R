library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(stringr)
library(data.table)



### load the key table

key_table<-read.csv(file = "key_table.csv") # refering to the modified key_table.csv that is used for mutect2 run script


### load the merged csv table containing the overview of all inserted mutations
sample_tables<-read.csv("merged_inserted_mutations_overview.csv", sep = ",") %>% select("csv_table", "position", "Variant", "occurrence") 
                        
colnames(sample_tables)<-c("sample_table", "position", "variant", "occurrence")

sample_tables$sample_table<-gsub("_config","",sample_tables$sample_table)

# add column for sample_tables identifier, fix the identifier to be identical to the one in the sample table
key_table$sample_table_id<-sub("-S[0-9]+","", key_table$sample_id)
key_table$sample_table_id<-gsub("-","_", key_table$sample_table_id)
key_table$sample_table_id <- gsub("(\\d)e_(\\d)", "\\1e-\\2", key_table$sample_table_id, perl=TRUE)
key_table$sample_table_id<-gsub("_tech_replicate_[0-9]","",key_table$sample_table_id)
key_table$sample_table_id<-paste0(key_table$sample_table_id, "_sample_table")


# filter sample_table table
sample_tables<-sample_tables %>% filter(sample_table %in% unique(key_table$sample_table_id))

## making a list of vcf_csv tables

normal_cancer <- key_table$normal_cancer

# Split the entries into numeric prefix and type
split_entries <- strsplit(normal_cancer, "_")
df <- data.frame(
  number_reads = sapply(split_entries, `[`, 1),
  number_copies = sapply(split_entries, `[`, 2),
  type=sapply(split_entries,`[`, 3),
  replicates=sapply(split_entries,`[`, 4),
  stringsAsFactors = FALSE
)

df[is.na(df)] <- 1

create_names <- function(df) {
  result <- c()
  
  # Unique combinations of number_reads and number_copies
  groups <- df %>%
    group_by(number_reads, number_copies) %>%
    group_split()
  
  for (group in groups) {
    # Extract the normal type row
    normal_row <- group %>% filter(type == "normal")
    
    # Combine normal with other types
    for (i in 1:nrow(group)) {
      if (group$type[i] != "normal") {
        name <- paste0(normal_row$number_reads, "_", normal_row$number_copies, "_normal_", 
                       group$number_reads[i], "_", group$number_copies[i], "_", 
                       group$type[i], "_", group$replicates[i])
        result <- c(result, name)
      }
    }
  }
  return(result)
}

# Create the names
final_results <- create_names(df)

# Print results
print(final_results)
final_results<-unique(final_results)
final_results<-paste0(final_results, ".csv")

key_table_short<-key_table %>% select(sample_index,read_number,copy_number,mutation_rate,normal_cancer,sample_table_id) %>% filter(!grepl("_normal",normal_cancer))


names_map <- setNames(final_results, sapply(final_results, function(x) {
  gsub(".csv", "", x)
}))

# Print the mapping to ensure it's correct
print(names_map)

names_map <- setNames(names_map, sub(".*normal_", "", sub(".csv$", "", names(names_map))))

# Initialize the new column with NA values
key_table_short$vcf_csv <- NA

# Iterate through each row of the DataFrame to match and assign the corresponding full filename
for (i in seq_along(key_table_short$normal_cancer)) {
  key_table_short$vcf_csv[i] <- names_map[key_table_short$normal_cancer[i]]
}

key_table_short$TP<-NA
key_table_short$FP<-NA
key_table_short$FN<-NA
key_table_short$precision<-NA
key_table_short$recall<-NA
key_table_short$F1<-NA



# Loop for TP, FP, FN, and precision and recall calculations.


for (n in 1:nrow(key_table_short)) {
  print(key_table_short$vcf_csv[n])
  vcf_table<-read.csv(paste0("r_processing/csv_tables/", key_table_short$vcf_csv[n]), sep = "\t")
  vcf_table<-vcf_table %>% select(2,4,5,7,9)
  vcf_table<-vcf_table %>% filter(FILTER == "PASS")
  
  sample_tables_sub<-sample_tables %>% filter(sample_table == key_table_short$sample_table_id[n])
  
  vcf_table_merged<-left_join(vcf_table,sample_tables_sub %>% select(position,variant,occurrence), by=c("POS" ="position", "ALT" = "variant"))
  
  TP=sum(!is.na(vcf_table_merged$occurrence))
  FP=sum(is.na(vcf_table_merged$occurrence))
  FN=nrow(sample_tables_sub)-TP
  
  precision=TP/(TP+FP)
  recall=TP/(TP+FN)
  
  F1=2*(precision*recall)/(precision+recall)
  
  key_table_short$TP[n]<-TP
  key_table_short$FP[n]<-FP
  key_table_short$FN[n]<-FN
  key_table_short$precision[n]<-precision
  key_table_short$recall[n]<-recall
  key_table_short$F1[n]<-F1
  
}

all_mut_info<-key_table_short %>% select(-sample_table_id, -vcf_csv)

all_mut_info <- all_mut_info %>%
  mutate(replicate = str_extract(normal_cancer, "(?<=_)[^_]+$"))

all_mut_info<-all_mut_info %>% select(sample_index, read_number,copy_number, mutation_rate, replicate, TP,FP,FN,precision, recall,F1)
all_mut_info$sample_index<-1:nrow(all_mut_info)

library(dplyr)
library(ggplot2)

all_mut_info[is.na(all_mut_info)]<-0

# Calculate summary statistics
summary_stats <- all_mut_info %>%
  group_by(read_number, copy_number, mutation_rate) %>%
  summarise(mean_precision = mean(precision),
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


summary_stats <- summary_stats %>%
  mutate(across(.cols = everything(), ~replace_na(., 0)))


max_precision <- max(summary_stats$mean_precision, na.rm = TRUE)
max_recall <- max(summary_stats$mean_recall, na.rm = TRUE)
desired_max_secondary_axis <- 1

# New scale factor based on desired max for the secondary axis
scale_factor <- max_precision / desired_max_secondary_axis


# Define your custom colors
custom_colors <- c("1500" = "#fdc827", "2500" = "#c23c81", "3500" = "#240691")

# Custom labels for x-axis
mutation_rates <- levels(factor(summary_stats$mutation_rate))
formatted_labels <- sapply(mutation_rates, function(x) {
  base <- 10
  exponent <- format(log10(as.numeric(x)), digits = 2)
  parse(text = paste0("paste('10^', '", exponent, "')"))
})




# Create a new column that combines read_number and copy_number for faceting
summary_stats$facet_label <- paste(summary_stats$read_number, summary_stats$copy_number, sep="_")

# Plotting
p <- ggplot(summary_stats, aes(x=factor(mutation_rate), group=factor(copy_number))) +
  geom_line(aes(y=mean_precision, color=factor(copy_number)), size=0.5, position=position_dodge(0.5)) +
  geom_point(aes(y=mean_precision), position=position_dodge(0.5), size=2, shape=21, fill="white") +
  geom_errorbar(aes(ymin=mean_precision-se_precision, ymax=mean_precision+se_precision),
                width=.1, position=position_dodge(0.5)) +
  geom_line(aes(y=mean_recall * scale_factor, color=factor(copy_number)), linetype="dashed", size=0.5, position=position_dodge(0.5)) +
  geom_point(aes(y=mean_recall * scale_factor), position=position_dodge(0.5), size=2, shape=21, fill="white", stroke=0.5) +
  geom_errorbar(aes(ymin=(mean_recall-se_recall)*scale_factor, ymax=(mean_recall+se_recall)*scale_factor),
                width=.1, position=position_dodge(0.5)) +
  facet_wrap(~facet_label, scales="free_x", nrow=1) +
  scale_color_manual(values=custom_colors) +
  theme_minimal() +
  labs(x="", y="") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none",
        panel.grid.major.y = element_blank(),  
        panel.grid.minor.y = element_blank())



# Save plot and summary table

ggsave("plot_precision_vs_recall_all_mut.jpg",p, height = 2, width = 7, units = "in")
write.csv(summary_stats, "summary_TP_FP_FN_rec_prec_run12.csv", row.names = F)
