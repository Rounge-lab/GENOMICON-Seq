library(Biostrings)
library(DBI)
library(RSQLite)
library(dplyr)
library(parallel)

###############################################################################################################
#                                                                                                             #
#      0.   W O R K I N G  D I R E C T O R Y  & S E T T I N G  A  S E E D  &  N U M B E R  O F  C O R E S     #
#                                                                                                             #
###############################################################################################################
setwd("./")

## Define the ouput directory
if (any(grepl("--output", commandArgs()))) { 
  output_arg <- which(grepl("--output", commandArgs())) 
  if (length(output_arg) > 0 && output_arg < length(commandArgs())) { 
    output_folder <- commandArgs()[output_arg + 1] 
    print(paste("Using", output_folder, "as the output folder"))
  } 
} 

# Create the output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

sql_database_folder<-"SQL_database"
if (!dir.exists(sql_database_folder)) {
  dir.create(sql_database_folder)
}



# argument --cores {integer} defines number of cores, if not given process is run on 1 core sequentially
if (any(grepl("--cores", commandArgs()))) {
  core_arg<-which(grepl("--cores", commandArgs()))
  if (length(core_arg) > 0 && core_arg < length(commandArgs())) {
    num_cores<-as.integer(commandArgs()[core_arg + 1])
    print(paste("Using", num_cores, "cores"))
  }
} else {
  num_cores <- 1
  print("Number of cores not provided, using deafult - 1 core. ALl processes will be run sequentially")
}

####################################################
#                                                  #
#      I.   D E F I N I N G   A R G U M E N T S    #
#                                                  #
####################################################

#Flag to control the flow
proceed_with_processing <- TRUE 

# Load the table with info about fasta file
if (any(grepl("^--fasta_lengths_info$", commandArgs()))) {
  # If fasta_sequence argument is provided, extract the file name
  fasta_arg_index <- which(grepl("^--fasta_lengths_info$", commandArgs(trailingOnly = TRUE)))
  
  if (fasta_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    fasta_content_file <- commandArgs(trailingOnly = TRUE)[fasta_arg_index + 1]
    fasta_content_path <- file.path(fasta_content_file)
    # Check if the file exists
    if (file.exists(fasta_content_path)) {
      # Read the file
      fasta_content <- read.csv(fasta_content_path)
      print("Fasta info loaded")
    } else {
      stop("The specified file does not exist.")
    }
  }
} else {
  stop("The fasta_lengths.cvs does not exist, please check the input directory")
}

# reading the path to fasta sequence
if (any(grepl("^--fasta_sequence$", commandArgs()))) {
  # If fasta_sequence argument is provided, extract the file name
  fasta_arg_index <- which(grepl("^--fasta_sequence$", commandArgs(trailingOnly = TRUE)))
  
  if (fasta_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    fasta_file <- commandArgs(trailingOnly = TRUE)[fasta_arg_index + 1]
    
    # Check if the file exists
    if (file.exists(fasta_file)) {
      print("fasta file is ready")
    } else {
      stop("The specified fasta file does not exist.")
    }
  }
} else {
  stop("Fasta sequence has to be provided, use --fasta_sequence /path/to/fasta.fa")
}

# Load the fasta sequence and check for existing SQLite files

  
    
if (exists("fasta_content")) {
  fasta_names<-fasta_content$fasta_names
  
  # Check for existing SQLite files
  existing_files <- sapply(fasta_names, function(name) {
    file.exists(file.path(sql_database_folder, paste0(name, ".sqlite")))
  })
  
  if (all(existing_files)) {
    print("All SQLite files already exist. Skipping processing.")
    proceed_with_processing <- FALSE #flag changes to FALSE
  } else {
    missing_sql_tables <- fasta_names[!existing_files]
    print(paste("Processing will continue for", missing_sql_tables, "chromosomes without existing SQLite files."))
    fasta_sequences<-readDNAStringSet(fasta_file)
    fasta_sequences<-fasta_sequences[missing_sql_tables]
  }
  
} else {
  stop("The specified fasta file does not exist.")
}



################################
#                              #
#      II F U N C T I O N S    #
#                              #
################################



process_chromosome <- function(chromosome, name) {
  sequence <- as.character(chromosome)
  nucleotides <- unlist(strsplit(sequence, ""))
  comp_nucleotides <- rev(chartr("ATCG", "TAGC", nucleotides))
  
  df <- data.frame(nucleotide = nucleotides,
                   comp_nucleotide = comp_nucleotides)
  
  # Tri and Penta context (Lead Strand)
  df$tri_context <- paste0(lag(df$nucleotide, 1, default = "-"), df$nucleotide, lead(df$nucleotide, 1, default = "-"))
  df$penta_context <- paste0(lag(df$nucleotide, 2, default = "-"), lag(df$nucleotide, 1, default = "-"), df$nucleotide, lead(df$nucleotide, 1, default = "-"), lead(df$nucleotide, 2, default = "-"))
  
  # Tri and Penta context (Reverse Strand)
  df$tri_context_rev_comp <- rev(paste0(lag(df$comp_nucleotide, 1, default = "-"), df$comp_nucleotide, lead(df$comp_nucleotide, 1, default = "-")))
  df$penta_context_rev_comp <- rev(paste0(lag(df$comp_nucleotide, 2, default = "-"), lag(df$comp_nucleotide, 1, default = "-"), df$comp_nucleotide, lead(df$comp_nucleotide, 1, default = "-"), lead(df$comp_nucleotide, 2, default = "-")))
  
  true_complement<-rev(comp_nucleotides)
  df$comp_nucleotide<-true_complement
  
  
  con <- dbConnect(SQLite(), file.path(sql_database_folder, paste0(name, ".sqlite")))
  dbWriteTable(con, name, df)
  dbDisconnect(con)
  cat("Processed chromosome:", name, "\n")
  return(name)
}


#############################################################################################################
#                                                                                                           #
#      III  R U N N I N G   T H E   F U N C T I O N    A N D   G E T T I N G   T H E   D A T A B A S E     #
#                                                                                                           #
#############################################################################################################


# Use mclapply to run the function in parallel for each chromosome if the files are not premade
if(proceed_with_processing) {
  results <- mclapply(1:length(fasta_sequences), function(i) {
    process_chromosome(fasta_sequences[[i]], names(fasta_sequences)[i])
  }, mc.cores = num_cores)
  
  print("Creating of SQL databases finished")
  
} else {
  print("No processing needed")
}


