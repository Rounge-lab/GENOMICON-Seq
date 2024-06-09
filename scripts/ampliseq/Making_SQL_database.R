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

# Load the fasta sequence and check for existing SQLite files
if (any(grepl("^--fasta_sequence$", commandArgs()))) {
  fasta_arg_index <- which(grepl("^--fasta_sequence$", commandArgs(trailingOnly = TRUE)))
  
  if (fasta_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    fasta_file <- commandArgs(trailingOnly = TRUE)[fasta_arg_index + 1]
    
    if (file.exists(fasta_file)) {
      fasta_content <- readDNAStringSet(fasta_file)
      print("fasta sequence loaded")
      
      # Check for existing SQLite files
      existing_files <- sapply(names(fasta_content), function(name) {
        file.exists(file.path(sql_database_folder, paste0(name, ".sqlite")))
      })
      
      if (all(existing_files)) {
        print("All SQLite files already exist. Skipping processing.")
        proceed_with_processing <- FALSE #flag changes to FALSE
      } else {
        fasta_content <- fasta_content[!existing_files]
        print(paste("Processing will continue for", length(fasta_content), "chromosomes without existing SQLite files."))
      }
      
    } else {
      stop("The specified fasta file does not exist.")
    }
  }
} else {
  stop("Fasta sequence has to be provided, use --fasta_sequence /path/to/fasta.fa")
}


## Is genome circular or not, by default is circular
if (any(grepl("--circular", commandArgs()))) {
  circular_pos <- which(grepl("--circular", commandArgs()))
  if (length(circular_pos) > 0 && circular_pos < length(commandArgs())) {
    circular_option <- commandArgs()[circular_pos + 1]
    if (circular_option == "Y") {
      print("Genome defined as circular")
    } else { 
      print("Genome defined as linear")}
  } 
} else{
  circular_option<- "Y"
  print("Genome not specified by user, by default is circular")
}


################################
#                              #
#      II F U N C T I O N S    #
#                              #
################################



process_chromosome <- function(chromosome, name, circular_option) {
  sequence <- as.character(chromosome)
  nucleotides <- unlist(strsplit(sequence, ""))
  comp_nucleotides <- rev(chartr("ATCG", "TAGC", nucleotides))
  
  if (circular_option == "Y") {
    # Extend for circular genome
    extended_nucleotides <- c(tail(nucleotides, 2), nucleotides, head(nucleotides, 2))
    extended_comp_nucleotides <- c(tail(comp_nucleotides, 2), comp_nucleotides, head(comp_nucleotides, 2))
  } else {
    extended_nucleotides <- nucleotides
    extended_comp_nucleotides <- comp_nucleotides
  }
  
  df <- data.frame(nucleotide = extended_nucleotides,
                   comp_nucleotide = extended_comp_nucleotides)
  
  # Tri and Penta context (Lead Strand)
  df$tri_context <- paste0(lag(df$nucleotide, 1, default = "-"), df$nucleotide, lead(df$nucleotide, 1, default = "-"))
  df$penta_context <- paste0(lag(df$nucleotide, 2, default = "-"), lag(df$nucleotide, 1, default = "-"), df$nucleotide, lead(df$nucleotide, 1, default = "-"), lead(df$nucleotide, 2, default = "-"))
  
  # Tri and Penta context (Reverse Strand)
  df$tri_context_rev_comp <- rev(paste0(lag(df$comp_nucleotide, 1, default = "-"), df$comp_nucleotide, lead(df$comp_nucleotide, 1, default = "-")))
  df$penta_context_rev_comp <- rev(paste0(lag(df$comp_nucleotide, 2, default = "-"), lag(df$comp_nucleotide, 1, default = "-"), df$comp_nucleotide, lead(df$comp_nucleotide, 1, default = "-"), lead(df$comp_nucleotide, 2, default = "-")))
  
  true_complement<-rev(extended_comp_nucleotides)
  df$comp_nucleotide<-true_complement
  
  if (circular_option == "Y") {
    # Remove the added nucleotides after calculating context
    df <- df[3:(nrow(df)-2), ]
  }
  
  
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
if( proceed_with_processing) {
  results <- mclapply(1:length(fasta_content), function(i) {
    process_chromosome(fasta_content[[i]], names(fasta_content)[i],circular_option)
  }, mc.cores = num_cores)
  
  print("Creating of SQL databases finished")
  
} else {
  print("No processing needed")
}


