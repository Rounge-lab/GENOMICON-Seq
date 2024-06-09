args <- commandArgs(trailingOnly = TRUE)
seed_value <- as.numeric(args[1])
num_seeds <- as.integer(args[2])

set.seed(seed_value)
seeds <- sample(1:1e8, num_seeds, replace = FALSE)

cat(seeds, sep="\n")

