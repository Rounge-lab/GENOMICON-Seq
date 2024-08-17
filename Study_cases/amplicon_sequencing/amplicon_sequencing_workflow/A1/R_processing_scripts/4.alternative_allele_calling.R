library("ggplot2")
library("ggpubr")
library("purrr")
library("dplyr")
library("tidyr")
#library("Hmisc")
library("tidyverse")
library(readxl)


#define functions
# Function for finding column name for specific value
maxN <- function(x, reachback = 0){
  # reachback = 0 is maximum, 1 is second to last, 2 is third to last and so on
  len <- length(x)
  if(reachback > len){
    error('You can not overreach the number of variables.')
  }
  names(sort(x, decreasing = TRUE)[1 + reachback])
}

# Function for finding the second largest value
maxN2 <- function(x, N = 2){
  len <- length(x)
  if(N > len){
    warning('N greater than length(x).  Setting N = length(x)')
    N <- length(x)
  }
  sort(x,partial = len - N+1)[len - N+1]
}


#coverage table
coverage_path <- c("A1_coverage_no_overhangs.csv")

coverage <- read.csv(paste(coverage_path), sep="\t")

# Filter out bases with QV <x
cutoff_qual<-30
coverage$A[coverage$qA < cutoff_qual] <- 0
coverage$T[coverage$qT < cutoff_qual] <- 0
coverage$C[coverage$qC < cutoff_qual] <- 0
coverage$G[coverage$qG < cutoff_qual] <- 0

# Separate F/R
Fvariants <- coverage[grep("-F", coverage$sample_id),]
Rvariants <- coverage[grep("-R", coverage$sample_id),]
Fvariants$sample_id<-gsub("-F","",Fvariants$sample_id)
Rvariants$sample_id<-gsub("-R","",Rvariants$sample_id)

#combine F and R
all_coverage=rbind(Fvariants, Rvariants) %>%
  group_by(sample_id,chr,position, reference) %>% mutate(CumCoverage=sum(coverage),
                                                          cumA = sum(A),
                                                          cumG =sum(G),
                                                          cumC = sum(C),
                                                          cumT = sum (T)) %>%
  dplyr::select(sample_id, chr, position, cumA, cumG, cumC, cumT,CumCoverage) %>%
  distinct()
#rename colnames
colnames(all_coverage) <- c("reference", "sample_id2", "chr", "position", "A", "G", "C", "T", "CumCoverage")



# call variants with functions
all_coverage[, "major_base"] <- apply(all_coverage[5:8], 1, function (x) maxN(x, reachback = 0))
all_coverage[, "major_base_count"] <- apply(all_coverage[, 5:8],1, max)
all_coverage[, "minor_base"] <- apply(all_coverage[5:8], 1, function (x) maxN(x, reachback = 1))
all_coverage[, "minor_base_count"] <- apply(all_coverage[, 5:8],1, function (x) maxN2(x))
all_coverage$coverage <- apply((all_coverage[5:8]), 1 ,sum)
all_coverage$minor_base_freq <- all_coverage$minor_base_count / all_coverage$coverage
all_coverage$minor_base_freq[is.nan(all_coverage$minor_base_freq)] <- 0
all_coverage$minor_base[all_coverage$minor_base_count == 0] <- "N"

all_coverage$reference <- toupper(all_coverage$reference)

#### Filter all variants whenever alternative nucleotide variant is >0

all_variants<-all_coverage %>% filter(minor_base_freq > 0)
write.csv(all_variants, "A1_all_variants.csv", row.names = FALSE)

