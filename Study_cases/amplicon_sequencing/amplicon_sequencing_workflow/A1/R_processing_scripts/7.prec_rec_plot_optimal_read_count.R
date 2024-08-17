library(ggplot2)
library(readr)
library(dplyr)

# Load the data - all mutations
data<-read.csv("A1_all_mutations_FP_TP_FN_summary_per_sample.csv")
data<-data %>% filter(optimal == "TRUE")
data<-data %>% select(original_sample_index,mutation_rate,polymerase_error_rate,mean_precision,mean_recall)

data$mutation_rate <- as.factor(data$mutation_rate)
data$polymerase_error_rate <- as.factor(data$polymerase_error_rate)
data$mean_precision <- as.numeric(data$mean_precision)
data$mean_recall <- as.numeric(data$mean_recall)

data$mutation_rate <- factor(data$mutation_rate, levels = c(5e-6, 5e-5, 5e-4))
data$polymerase_error_rate <- as.factor(data$polymerase_error_rate)

plot <- ggplot(data, aes(x = polymerase_error_rate, group = 1)) +
  geom_point(aes(y = mean_precision, color = "mean_prec"), size = 3) + # Change size of the dots
  geom_point(aes(y = mean_recall, color = "mean_rec"), size = 3) + # Change size of the dots
  geom_line(aes(y = mean_precision, color = "mean_prec"), linetype = "dashed") + # Make lines dashed
  geom_line(aes(y = mean_recall, color = "mean_rec"), linetype = "dashed") + # Make lines dashed
  facet_wrap(~ mutation_rate, scales = "free_x", strip.position = "top") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("mean_prec" = "#4F1B77", "mean_rec" = "#DF6665")) + # Customize colors here
  labs(x = "Pol Rate", y = "Mean precision/recall", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "none", # Remove legend
        panel.grid.major.y = element_blank(), # Remove horizontal major grid lines
        panel.grid.minor.y = element_blank()) # Remove horizontal minor grid lines



ggsave("A1_summarized_prec_rec_optimal_read_count.jpg", plot = plot, width = 8, height = 2.5)


# Load the data - APOBEC3
dataA3<-read.csv("A1_APOBEC_mutations_FP_TP_FN_summary_per_sample.csv")
dataA3<-dataA3 %>% filter(optimal == "TRUE")
dataA3<-dataA3 %>% select(original_sample_index,mutation_rate,polymerase_error_rate,mean_precision,mean_recall)

dataA3$mutation_rate <- as.factor(dataA3$mutation_rate)
dataA3$polymerase_error_rate <- as.factor(dataA3$polymerase_error_rate)
dataA3$mean_precision <- as.numeric(dataA3$mean_precision)
dataA3$mean_recall <- as.numeric(dataA3$mean_recall)

dataA3$mutation_rate <- factor(dataA3$mutation_rate, levels = c(5e-6, 5e-5, 5e-4))
dataA3$polymerase_error_rate <- as.factor(dataA3$polymerase_error_rate)

plot2 <- ggplot(dataA3, aes(x = polymerase_error_rate, group = 1)) +
  geom_point(aes(y = mean_precision, color = "mean_prec"), size = 3) + # Change size of the dots
  geom_point(aes(y = mean_recall, color = "mean_rec"), size = 3) + # Change size of the dots
  geom_line(aes(y = mean_precision, color = "mean_prec"), linetype = "dashed") + # Make lines dashed
  geom_line(aes(y = mean_recall, color = "mean_rec"), linetype = "dashed") + # Make lines dashed
  facet_wrap(~ mutation_rate, scales = "free_x", strip.position = "top") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("mean_prec" = "#4F1B77", "mean_rec" = "#DF6665")) + # Customize colors here
  labs(x = "Pol Rate", y = "Mean precision/recall", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "none", # Remove legend
        panel.grid.major.y = element_blank(), # Remove horizontal major grid lines
        panel.grid.minor.y = element_blank()) # Remove horizontal minor grid lines



ggsave("A1_summarized_prec_rec_optimal_read_count_APOBEC3.jpg", plot = plot2, width = 8, height = 2.5)