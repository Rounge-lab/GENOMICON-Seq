library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(stringr)
library(data.table)

#### load files
original_mutations<-read.csv("merged_inserted_mutations_overview.csv", sep = ",") %>% select("csv_table", "position", "Variant", "occurrence")
colnames(sample_tables)<-c("sample_table", "position", "variant", "occurrence")

detected_muts<-read.csv("35M_2500_normal_35M_2500_0.4_sbs24.csv", sep = "\t") # vcf files converetd into csv table, here is an example of sample with sbs24 signature
detected_muts<-detected_muts %>% filter(FILTER == "PASS")
detected_muts<-detected_muts %>% select(POS,REF,ALT)

# extracted positions

row_ids<-detected_muts$POS

# Query the SQL-database
# Use the SQL-database that is part of the genomicon-seq
library(DBI)
library(RSQLite)

conn <- dbConnect(RSQLite::SQLite(), dbname = "SQL_database/chr1.sqlite")

ids_string <- paste(row_ids, collapse = ", ")

query <- sprintf("SELECT ROWID AS positions, * FROM chr1 WHERE ROWID IN (%s)", ids_string)
result <- dbGetQuery(conn, query)
dbDisconnect(conn)

### now fix the table so that there are 16 possible tricontexts

passed_mutations<-result %>% select(positions,nucleotide,comp_nucleotide,tri_context,tri_context_rev_comp)

passed_mutations_merged<-left_join(passed_mutations, detected_muts, by=c("positions" = "POS", "nucleotide" = "REF"))

passed_mutations_merged$ALT_rev_comp <- ifelse(passed_mutations_merged$ALT == "A", "T",
                          ifelse(passed_mutations_merged$ALT == "T", "A",
                                 ifelse(passed_mutations_merged$ALT == "C", "G", "C")))

all_contexts<-c("ACA", "ACT", "ACC", "ACG", "CCA", "CCT", "CCC", "CCG", "GCA", "GCT", "GCC", "GCG", "TCA", "TCT", "TCC", "TCG","ATA", "ATT", "ATC", "ATG", "CTA", "CTT", "CTC", "CTG", "GTA", "GTT", "GTC", "GTG", "TTA", "TTT", "TTC", "TTG")


passed_mutations_merged$context <- apply(passed_mutations_merged, 1, function(x) {
  context_match <- all_contexts[all_contexts %in% c(x["tri_context"], x["tri_context_rev_comp"])]
  if (length(context_match) > 0) {
    return(context_match[1])  # Returning the first match if there are multiple
  } else {
    return(NA)  # Return NA if no match found
  }
})

passed_mutations_merged$nuc_context<-ifelse(passed_mutations_merged$context == passed_mutations_merged$tri_context, passed_mutations_merged$nucleotide,passed_mutations_merged$comp_nucleotide)

passed_mutations_merged$substituted_nuc<-ifelse(passed_mutations_merged$nuc_context == passed_mutations_merged$nucleotide, passed_mutations_merged$ALT, passed_mutations_merged$ALT_rev_comp)

### filter out those with more then 1 char in substituted nuc
passed_mutations_merged<-passed_mutations_merged %>% filter(nchar(substituted_nuc) == 1)

passed_mutations_merged$sustitution_type<-paste0(passed_mutations_merged$nuc_context,">",passed_mutations_merged$substituted_nuc)

passed_mutations_short<-passed_mutations_merged %>% select(positions,context,sustitution_type)

## Now the orignal mutations

row_ids<-original_mutations$position

# Query the SQL database
library(DBI)
library(RSQLite)

conn <- dbConnect(RSQLite::SQLite(), dbname = "SQL_database/chr1.sqlite")

ids_string <- paste(row_ids, collapse = ", ")

query <- sprintf("SELECT ROWID AS positions, * FROM chr1 WHERE ROWID IN (%s)", ids_string)
result <- dbGetQuery(conn, query)
dbDisconnect(conn)


original_mutations_all<-result %>% select(positions,nucleotide,comp_nucleotide,tri_context,tri_context_rev_comp)

original_mutations_all_merged<-left_join(original_mutations_all, original_mutations, by=c("positions" = "position"))

original_mutations_all_merged$variant_rev_comp <- ifelse(original_mutations_all_merged$variant == "A", "T",
                                               ifelse(original_mutations_all_merged$variant == "T", "A",
                                                      ifelse(original_mutations_all_merged$variant == "C", "G", "C")))


original_mutations_all_merged$context <- apply(original_mutations_all_merged, 1, function(x) {
  context_match <- all_contexts[all_contexts %in% c(x["tri_context"], x["tri_context_rev_comp"])]
  if (length(context_match) > 0) {
    return(context_match[1])  # Returning the first match if there are multiple
  } else {
    return(NA)  # Return NA if no match found
  }
})

original_mutations_all_merged$nuc_context<-ifelse(original_mutations_all_merged$context == original_mutations_all_merged$tri_context, original_mutations_all_merged$nucleotide,original_mutations_all_merged$comp_nucleotide)

original_mutations_all_merged$substituted_nuc<-ifelse(original_mutations_all_merged$nuc_context == original_mutations_all_merged$nucleotide, original_mutations_all_merged$variant, original_mutations_all_merged$variant_rev_comp)


original_mutations_all_merged$sustitution_type<-paste0(original_mutations_all_merged$nuc_context,">",original_mutations_all_merged$substituted_nuc)


original_mutations_short<-original_mutations_all_merged %>% select(positions,context,sustitution_type, occurance)



#### fidning TP FN FP


# Find True Positives (TPs)
tp <- merge(passed_mutations_short, original_mutations_short, by = c("positions", "context", "sustitution_type"))

# Find False Positives (FPs) - in table1 but not in tp_table
fp <- setdiff(passed_mutations_short, tp_table[, names(passed_mutations_short)])

# Find False Negatives (FNs) - in table2 but not in tp_table
fn <- setdiff(original_mutations_short, tp_table[, names(original_mutations_short)])

# Print the results
print("True Positives (TPs):")
print(nrow(tp))
print("False Positives (FPs):")
print(nrow(fp))
print("False Negatives (FNs):")
print(nrow(fn))


## plot passed_mutations_short

summary_passed<-passed_mutations_short %>%
  group_by(context,sustitution_type) %>%
  summarise(count_positions = n())


## make template 

substitutions<-c("C>A" , "C>G", "C>T", "T>A", "T>C", "T>G")

# Filter contexts for C or T in the middle position
contexts_with_c <- all_contexts[grep("^.C.$", all_contexts)]
contexts_with_t <- all_contexts[grep("^.T.$", all_contexts)]

# Filter substitutions starting with C or T
subs_with_c <- substitutions[grep("^C", substitutions)]
subs_with_t <- substitutions[grep("^T", substitutions)]

# Expand each context with the corresponding substitutions
expanded_with_c <- expand.grid(context = contexts_with_c, substitution = subs_with_c)
expanded_with_t <- expand.grid(context = contexts_with_t, substitution = subs_with_t)

# Combine both expansions
final_table <- rbind(expanded_with_c, expanded_with_t)

summary_passed_all<-left_join(final_table,summary_passed, by=c("context" = "context", "substitution" = "sustitution_type"))
summary_passed_all[is.na(summary_passed_all)] <- 0

summary_passed_all$percentage<- summary_passed_all$count_positions / sum(summary_passed_all$count_positions) * 100

# Define colors for each substitution type
substitution_colors <- c(
  "C>A" = "#1f78b4", # Blue
  "C>G" = "#000000", # Black
  "C>T" = "#e31a1c", # Red
  "T>A" = "#C9C8C9", # Light blue
  "T>C" = "#b2df8a", # Green
  "T>G" = "#fb9a99"  # Pink
)

# Plot
plot_after_pos_counts<-ggplot(summary_passed_all, aes(x = context, y = percentage, fill = substitution)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ substitution, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = substitution_colors) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
  ) +
  labs(
    y = "Percentage of Single Base Substitutions",
  ) + theme(panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank())+
  scale_y_continuous(limits = c(0, 15))

plot_after_pos_counts
ggsave("sbs24_post_sequencinng_counts.jpg", plot_after_pos_counts,width = 11, height = 2, units = "in")


### summary for pre sequencing mutation positions


summary_original<-original_mutations_short %>%
  group_by(context,sustitution_type) %>%
  summarise(count_positions = n())

summary_original_all<-left_join(final_table,summary_original, by=c("context" = "context", "substitution" = "sustitution_type"))
summary_original_all[is.na(summary_original_all)] <- 0
summary_original_all$percentage<-summary_original_all$count_positions / sum(summary_original_all$count_positions) * 100

plot_before_pos_counts<-ggplot(summary_original_all, aes(x = context, y = percentage, fill = substitution)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ substitution, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = substitution_colors) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
  ) +
  labs(
    y = "Percentage of Single Base Substitutions",
  )+ theme(panel.grid.major.x = element_blank(), 
           panel.grid.minor.x = element_blank(),
           panel.grid.minor.y = element_blank())+
  scale_y_continuous(limits = c(0, 15))

plot_before_pos_counts
ggsave("sbs24_pre_sequencinng_counts.jpg", plot_before_pos_counts,width = 11, height = 2, units = "in")

### plot the original sbs signtures - shorter csv table made from the orignal sbs-table downloaded from COSMIC

sbs_original<-read.csv(file = "sbs24_short.csv")
transformed_original<-sbs_original %>%
  mutate(
    Mutation = str_extract(X, "\\[.*?\\]"),         # Extract the mutation pattern including brackets
    Mutation = str_remove_all(Mutation, "\\[|\\]"), # Remove the brackets
    X = str_replace_all(X, "\\[.*?\\]", substr(Mutation, 1, 1)), # Replace the mutation pattern with the first character
    X = str_replace_all(X, ">$", "")  # Remove the 'greater than' sign at the end if exists
  )

transformed_original<-transformed_original %>% select(X,Mutation,SBS24_GRCh38)

transformed_original<-transformed_original %>% rename(context = X, substitution=Mutation, percentage=SBS24_GRCh38)
transformed_original$percentage<-transformed_original$percentage*100

# Plot
original_data<-ggplot(transformed_original, aes(x = context, y = percentage, fill = substitution)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ substitution, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = substitution_colors) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    panel.spacing.x = unit(0, "lines"),
  ) +
  labs(
    y = "Percentage of Single Base Substitutions",
  ) + theme(panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank())+
  scale_y_continuous(limits = c(0, 15))

original_data
ggsave("sbs24_original_data.jpg", original_data,width = 11, height = 2, units = "in")
