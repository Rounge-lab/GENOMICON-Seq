###############################################################################
#  Fragmentation of genomes into random fragments (linear or circular)       #
###############################################################################

library(Biostrings)
library(foreach)
library(iterators)
library(doParallel)
library(data.table)
library(Rcpp)

#########################################################################
# I.  PARSE ARGUMENTS                                                    #
#########################################################################

args <- commandArgs()

# --input_directory
input_folder <- if (any(grepl("--input_directory", args))) {
  args[which(grepl("--input_directory", args)) + 1]
} else stop("input_directory must be supplied")

# --output_directory
output_folder <- if (any(grepl("--output_directory", args))) {
  args[which(grepl("--output_directory", args)) + 1]
} else stop("output_directory must be supplied")

setwd(input_folder)

# --set_seed
if (any(grepl("--set_seed", args))) {
  seed_value <- as.integer(args[which(grepl("--set_seed", args)) + 1])
  if (is.na(seed_value) || seed_value < 0) stop("Invalid seed")
  set.seed(seed_value)
  message("Seed set to ", seed_value)
}

# --fragmentation_fraction
fragmentation_fraction <- if (any(grepl("--fragmentation_fraction", args))) {
  as.numeric(args[which(grepl("--fragmentation_fraction", args)) + 1])
} else 1
if (is.na(fragmentation_fraction) || fragmentation_fraction <= 0 || fragmentation_fraction > 1)
  stop("fragmentation_fraction must be in (0,1]")
message("Fragmentation fraction: ", fragmentation_fraction)

# --circular  (Y/N – default Y)
circular_option <- if (any(grepl("--circular", args))) {
  args[which(grepl("--circular", args)) + 1]
} else "Y"
circular_option <- toupper(circular_option)
if (!circular_option %in% c("Y","N"))
  stop("--circular must be Y or N")
message("Genome treated as ", ifelse(circular_option=="Y","circular","linear"))

# --fragment_length  (single int or min:max)
if (any(grepl("--fragment_length", args))) {
  fl_arg <- args[which(grepl("--fragment_length", args))+1]
  if (grepl(":", fl_arg)) {
    range_vals <- as.integer(strsplit(fl_arg, ":")[[1]])
    if (length(range_vals)!=2 || any(is.na(range_vals)) || any(range_vals<=0) ||
        range_vals[1]>range_vals[2]) stop("Invalid fragment_length range")
    fragment_length <- range_vals[1]:range_vals[2]
  } else {
    fl <- as.integer(fl_arg)
    if (is.na(fl) || fl<=0) stop("Invalid fragment_length")
    fragment_length <- fl
  }
} else fragment_length <- 250:450
message("Fragment length(s): ", if (length(fragment_length)==1) fragment_length
                               else paste0(min(fragment_length),":",max(fragment_length)))

# --cores
num_cores <- if (any(grepl("--cores", args))) {
  as.integer(args[which(grepl("--cores", args))+1])
} else 2
if (is.na(num_cores) || num_cores<=0) stop("Invalid --cores")
message("Using ", num_cores, " core(s)")

#########################################################################
# II.  FUNCTIONS                                                        #
#########################################################################

expand_row <- function(row, fasta_name_length) {
  genome_name_base <- row$genome_name
  copy_number      <- row$copy_number
  genome_len       <- fasta_name_length[[sub("_.*","", genome_name_base)]]
  data.table(genome_name = paste0(genome_name_base, "_c", seq_len(copy_number)),
             new_start   = sample.int(genome_len, copy_number, replace = TRUE))
}

expand_data_table <- function(dt, fasta_name_length, seed) {
  set.seed(seed)
  rbindlist(lapply(seq_len(nrow(dt)), function(i) expand_row(dt[i], fasta_name_length)))
}

apply_binomial <- function(dt, column_name, eff) {
  dt[[column_name]] <- rbinom(nrow(dt), dt[[column_name]], prob = eff)
  dt
}

#########################################################################
# III.  LOAD SAMPLE TABLES & FASTA LENGTHS                              #
#########################################################################

sample_table_path <- list.files(pattern = "_sample_table\\.csv$")
if (!length(sample_table_path)) stop("No *_sample_table.csv found")

fasta_name_table  <- fread("fasta_lengths.csv")
fasta_name_length <- split(fasta_name_table$fasta_lengths, fasta_name_table$fasta_names)

data_tables <- setNames(lapply(sample_table_path, fread),
                        sub("\\.csv$","", sample_table_path))

data_tables <- lapply(data_tables, apply_binomial,
                      column_name = "copy_number",
                      eff = fragmentation_fraction)

data_tables <- mclapply(data_tables, function(dt) dt[copy_number>0],
                        mc.cores = num_cores)

# per-worker seeds
worker_seeds <- sample.int(1e7, length(data_tables))

#########################################################################
# IV.  NEW STARTS (circular genomes)                                    #
#########################################################################

if (circular_option=="Y") {
  new_starts <- mclapply(seq_along(data_tables), function(i)
    expand_data_table(data_tables[[i]], fasta_name_length, worker_seeds[i]),
    mc.cores = num_cores)
  names(new_starts) <- names(data_tables)
} else {
  new_starts <- lapply(data_tables, function(dt) {
    dt2 <- dt[rep(seq_len(nrow(dt)), dt$copy_number)]
    dt2[, new_start := 1L]
    dt2[, genome_name := paste0(genome_name,"_c", seq_len(.N)), by = genome_name]
  })
}

#########################################################################
# V.  FRAGMENT COORDINATES VIA C++                                      #
#########################################################################

Rcpp::sourceCpp("/usr/src/app/scripts/ampliseq/fragmentation.cpp")

mclapply(seq_along(data_tables), function(i) {
  tbl_name   <- names(data_tables)[i]
  prefix     <- sub("_.*","", tbl_name)
  genome_len <- fasta_name_length[[prefix]][1]
  expand_table_cpp(data_tables[[i]], tbl_name, genome_len,
                   min(fragment_length), max(fragment_length),
                   worker_seeds[i],
                   circular_option=="Y",          #  <-- NEW FLAG
                   output_folder)
}, mc.cores = num_cores)

#########################################################################
# VI.  WRITE helper CSVs                                                #
#########################################################################

save_tables <- function(tables, suffix) {
  names(tables) <- sub("_sample_table", suffix, names(tables))
  for (nm in names(tables)) {
    fwrite(tables[[nm]], file.path(output_folder, paste0(nm, ".csv")))
  }
}

save_tables(data_tables, "_fragmentation_selected")
save_tables(new_starts,  "_new_start")

message("Fragmentation finished!")
