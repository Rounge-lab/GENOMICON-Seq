library(ggplot2)
library(readr)
library(dplyr)


data2<-read.csv("A2_all_mutations_FP_TP_FN_summary_per_sample.csv")
data2<-data2 %>% filter(optimal == "TRUE")
data2<-data2%>% select(original_sample_index,number_copies,fragment_fraction,mean_precision,mean_recall)
data2$number_copies<-as.character(data2$number_copies)

data2$number_copies<-factor(data2$number_copies, levels = c("10000","100000","200000"))
data2$fragment_fraction<-as.factor(data2$fragment_fraction)
data2$mean_precision<-as.numeric(data2$mean_precision)
data2$mean_recall<-as.numeric(data2$mean_recall)

# Prepare the plot with custom colors
plot <- ggplot(data2, aes(x = fragment_fraction, group = 1)) +
  geom_point(aes(y = mean_precision, color = "mean_prec"), size = 3) + # Change size of the dots
  geom_point(aes(y = mean_recall, color = "mean_rec"), size = 3) + # Change size of the dots
  geom_line(aes(y = mean_precision, color = "mean_prec"), linetype = "dashed") + # Make lines dashed
  geom_line(aes(y = mean_recall, color = "mean_rec"), linetype = "dashed") + # Make lines dashed
  facet_wrap(~ number_copies, scales = "free_x", strip.position = "top") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("mean_prec" = "#4F1B77", "mean_rec" = "#DF6665")) + # Customize colors here
  labs(x = "sample input", y = "Mean precision/recall", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "none", # Remove legend
        panel.grid.major.y = element_blank(), # Remove horizontal major grid lines
        panel.grid.minor.y = element_blank()) # Remove horizontal minor grid lines

# Print the plot
print(plot)

ggsave("A2_summarized_prec_rec_optimal_read_count.jpg", plot = plot, width = 6, height = 2.5)


# Load the data - APOBEC3

data2A3<-read.csv("A2_APOBEC_mutations_FP_TP_FN_summary_per_sample.csv")
data2A3<-data2A3 %>% filter(optimal == "TRUE")
data2A3<-data2A3%>% select(original_sample_index,number_copies,fragment_fraction,mean_precision,mean_recall)
data2A3$number_copies<-as.character(data2A3$number_copies)

data2A3$number_copies<-factor(data2A3$number_copies, levels = c("10000","100000","200000"))
data2A3$fragment_fraction<-as.factor(data2A3$fragment_fraction)
data2A3$mean_precision<-as.numeric(data2A3$mean_precision)
data2A3$mean_recall<-as.numeric(data2A3$mean_recall)

# Prepare the plot with custom colors
plot2 <- ggplot(data2A3, aes(x = fragment_fraction, group = 1)) +
  geom_point(aes(y = mean_precision, color = "mean_prec"), size = 3) + # Change size of the dots
  geom_point(aes(y = mean_recall, color = "mean_rec"), size = 3) + # Change size of the dots
  geom_line(aes(y = mean_precision, color = "mean_prec"), linetype = "dashed") + # Make lines dashed
  geom_line(aes(y = mean_recall, color = "mean_rec"), linetype = "dashed") + # Make lines dashed
  facet_wrap(~ number_copies, scales = "free_x", strip.position = "top") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("mean_prec" = "#4F1B77", "mean_rec" = "#DF6665")) + # Customize colors here
  labs(x = "sample input", y = "Mean precision/recall", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "none", # Remove legend
        panel.grid.major.y = element_blank(), # Remove horizontal major grid lines
        panel.grid.minor.y = element_blank()) # Remove horizontal minor grid lines

# Print the plot
print(plot2)

ggsave("A2_summarized_prec_rec_optimal_read_count_APOBEC3.jpg", plot = plot2, width = 6, height = 2.5)


