
### Combine read count and coverage information

rm(list=ls())

library("ggplot2")
library("reshape2")
library("splitstackshape")
library("tidyr")
library("dplyr")


### Upload files
#compiled stats
stats <- read.csv("A3_compiled_stats.csv", header = T)
#no overhangs coverage
coverage <- read.csv("A3_coverage_no_overhangs.csv",header = T, sep = "\t") 

### Create output file
stats$sample_id <- gsub("-F", "", stats$sample_id)
stats$sample_id <- gsub("-R", "", stats$sample_id)
stats$sample_id <- gsub("-S.*$", "", stats$sample_id)

stats_sum <- stats %>%
  group_by(sample_id, strain, reference) %>%
    summarise(raw_reads=sum(raw_reads),
              trimmed_reads=sum(trimmed_reads)
          )


coverage$chr <-  str_extract(coverage$chr, "HPV[0-9]+")


coverageF <- coverage[grep("-F", coverage$sample_id), ]
coverageR <- coverage[grep("-R", coverage$sample_id), ]

coverageF$sample_id<-gsub("-F","",coverageF$sample_id)
coverageR$sample_id<-gsub("-R","",coverageR$sample_id)


all_coverage=rbind(coverageF, coverageR) %>%
   group_by(sample_id, position) %>% mutate(coverage_sum=sum(coverage)) %>%
   dplyr::select(sample_id, chr, position, coverage_sum) %>%
   distinct()

write.csv(all_coverage, "A3_cum_cov_per_site.csv")

genomecov <- all_coverage %>%
  group_by(sample_id, chr) %>%
  summarise(genomelen = n(),
            mean(coverage_sum),
            median(coverage_sum),
            min(coverage_sum),
            max(coverage_sum),
            count10 = sum(coverage_sum > 9),           
            coverage10 = count10 / genomelen,
            count50 = sum(coverage_sum > 49),           
            coverage50 = count50 / genomelen,
            count100 = sum(coverage_sum > 99),           
            coverage100 = count100 / genomelen
          )

genomecov <- genomecov[, !(colnames(genomecov) %in% c("count10","count50","count100"))]
# genomecov$sample_id <- gsub("-S01","", genomecov$sample_id)

colnames(genomecov) <- c("sample_id","hpv","genome_length","mean_coverage","median_coverage","min_coverage","max_coverage",
                               "genome_covered_10x","genome_covered_50x","genome_covered_100x")

count_coverage <- merge(stats_sum, genomecov,by="sample_id",all=T)
count_coverage <- count_coverage[, !(colnames(count_coverage) %in% c("hpv","genome_length"))]

write.csv(file="A3_counts_coverage.csv", count_coverage, row.names = FALSE)


