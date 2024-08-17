library(ggplot2)
library(readr)
library(dplyr)

data3<-read.csv("A1_all_mutations_FP_TP_FN_summary_per_sample.csv")

data3<-data3 %>% filter(optimal == "TRUE")
data3$number_reads<-c("100k", "500k", "1.5M", "5M")
data3<-data3%>% select(original_sample_index,number_reads,mean_precision,mean_recall)
data3$number_reads<-factor(data3$number_reads, levels = c("100k","500k","1.5M","5M"))

data3$mean_precision<-as.numeric(data3$mean_precision)
data3$mean_recall<-as.numeric(data3$mean_recall)

plot3 <- ggplot(data3, aes(x = number_reads, group = 1)) +
  geom_point(aes(y = mean_precision, color = "mean_prec"), size = 3.6) + # Change size of the dots
  geom_point(aes(y = mean_recall, color = "mean_rec"), size = 3.6) + # Change size of the dots
  geom_line(aes(y = mean_precision, color = "mean_prec"), linetype = "dashed") + # Make lines dashed
  geom_line(aes(y = mean_recall, color = "mean_rec"), linetype = "dashed") + # Make lines dashed
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("mean_prec" = "#4F1B77", "mean_rec" = "#DF6665")) + # Customize colors here
  labs(x = "Pol Rate", y = "Values", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "none", # Remove legend
        panel.grid.major.y = element_blank(), # Remove horizontal major grid lines
        panel.grid.minor.y = element_blank()) # Remove horizontal minor grid lines

ggsave("A3_summarized_prec_rec_optimal_read_count.jpg", plot = plot3, width = 3, height = 2.5)


#APOBEC3

data3A3<-read.csv("A1_APOBEC_FP_TP_FN_summary_per_sample.csv")

data3A3<-data3A3 %>% filter(optimal == "TRUE")
data3A3$number_reads<-c("100k", "500k", "1.5M", "5M")
data3A3<-data3A3%>% select(original_sample_index,number_reads,mean_precision,mean_recall)
data3A3$number_reads<-factor(data3A3$number_reads, levels = c("100k","500k","1.5M","5M"))

data3A3$mean_precision<-as.numeric(data3A3$mean_precision)
data3A3$mean_recall<-as.numeric(data3A3$mean_recall)

plot3 <- ggplot(data3A3, aes(x = number_reads, group = 1)) +
  geom_point(aes(y = mean_precision, color = "mean_prec"), size = 3.6) + # Change size of the dots
  geom_point(aes(y = mean_recall, color = "mean_rec"), size = 3.6) + # Change size of the dots
  geom_line(aes(y = mean_precision, color = "mean_prec"), linetype = "dashed") + # Make lines dashed
  geom_line(aes(y = mean_recall, color = "mean_rec"), linetype = "dashed") + # Make lines dashed
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = c("mean_prec" = "#4F1B77", "mean_rec" = "#DF6665")) + # Customize colors here
  labs(x = "Pol Rate", y = "Values", color = "Legend") +
  theme_minimal() +
  theme(legend.position = "none", # Remove legend
        panel.grid.major.y = element_blank(), # Remove horizontal major grid lines
        panel.grid.minor.y = element_blank()) # Remove horizontal minor grid lines

ggsave("A3_summarized_prec_rec_optimal_read_count_APOBEC3.jpg", plot = plot3, width = 3, height = 2.5)
