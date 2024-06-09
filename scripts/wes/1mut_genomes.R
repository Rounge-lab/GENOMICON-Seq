library(dplyr)
library(stringr)
library("Biostrings")
library(stringi)
library(DBI)
library(RSQLite)
library(data.table)
library(foreach)
library(iterators)
library(doParallel)
library(tidyr)
###############################################################################################################
#                                                                                                             #
#      0.   W O R K I N G  D I R E C T O R Y  & S E T T I N G  A  S E E D  &  N U M B E R  O F  C O R E S     #
#                                                                                                             #
###############################################################################################################

if (any(grepl("--output", commandArgs()))) { 
  output_arg <- which(grepl("--output", commandArgs())) 
  if (length(output_arg) > 0 && output_arg < length(commandArgs())) { 
    output_folder <- commandArgs()[output_arg + 1] 
    print(paste("Using", output_folder, "as the output folder")) } } 


if (any(grepl("--set_seed", commandArgs()))) {
  seed_arg <- which(grepl("--set_seed", commandArgs()))
  if (length(seed_arg) > 0 && seed_arg < length(commandArgs())) {
    seed_value <- as.integer(commandArgs()[seed_arg + 1])
    if (!is.na(seed_value) && seed_value >= 0) {
      set.seed(seed_value)
      print(paste0("Seed set to ", seed_value))
    } else {
      stop("Invalid seed value")
    }
  }
} else {
  print("Seed not provided")
}


if (any(grepl("--cores", commandArgs()))) {
  core_arg <- which(grepl("--cores", commandArgs()))
  if (length(core_arg) > 0 && core_arg < length(commandArgs())) {
    num_cores <- as.integer(commandArgs()[core_arg + 1])
    
    if (is.na(num_cores) || num_cores <= 0 || num_cores != floor(num_cores)) {
      stop("Invalid number of cores. It must be a positive integer.")
    }
    
    available_cores <- parallel::detectCores()
    if (num_cores > available_cores) {
      warning(paste("Requested number of cores (", num_cores, ") exceeds the available cores (", available_cores, "). Proceed with caution."))
    }
    
    print(paste("Using", num_cores, "cores"))
  }
} else {
  num_cores <- 1
  print("Number of cores not provided, using default - 1 core. All processes will be run sequentially")
}


##################################################
#                                                #
#      I.   D E F I N I N G   S E Q U E N C E    #
#                                                #
##################################################

# Load the fasta sequence
if (any(grepl("^--fasta_lengths_info$", commandArgs()))) {
  fasta_arg_index <- which(grepl("^--fasta_lengths_info$", commandArgs(trailingOnly = TRUE)))
  
  if (fasta_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    fasta_content_file <- commandArgs(trailingOnly = TRUE)[fasta_arg_index + 1]
    fasta_content_path <- file.path(fasta_content_file)
    if (file.exists(fasta_content_path)) {
      fasta_content <- read.csv(fasta_content_path)
      print("Fasta info loaded")
    } else {
      stop("The specified file does not exist.")
    }
  }
} else {
  stop("The fasta_lengths.cvs does not exist, please check the input directory")
}

# Load the probe targeted region (exom bed)

if (any(grepl("^--target_exome_region$", commandArgs()))) {
  exome_arg_index <- which(grepl("^--target_exome_region$", commandArgs(trailingOnly = TRUE)))
  
  if (exome_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    bed_data_file <- commandArgs(trailingOnly = TRUE)[exome_arg_index + 1]
    bed_data_path <- file.path(bed_data_file)
    # Check if the file exists
    if (file.exists(bed_data_path)) {
      # Read the file
      bed_data <- read.table(bed_data_path, header=FALSE, sep="\t")
      print("Bed file loaded")
    } else {
      stop("The specified file does not exist.")
    }
  }
} else {
  stop("The bed does not exist, please check the input directory")
}


################################################################################################
#                                                                                              #
#      II. N U C L E O T I D E  C O N T E X T  &  M U T A T I O N  S I G N A T U R E           #
#                                                                                              #
################################################################################################
if (any(grepl("--mutational_context", commandArgs()))) {
  print("Mutation context is defined")
  
  mut_context_pos <- which(grepl("--mutational_context", commandArgs()))
  if (length(mut_context_pos) > 0 && mut_context_pos < length(commandArgs())) {
    mut_context <- unlist(strsplit(commandArgs()[mut_context_pos + 1], ","))
    
    # Check if all elements have equal length and are either 1, 3, or 5 characters long
    char_lengths <- nchar(mut_context)
    all_same_length <- all(char_lengths == char_lengths[1])
    valid_lengths <- all(char_lengths %in% c(1, 3, 5))
    
    if (all_same_length && valid_lengths) {
      print(mut_context)
    } else {
      stop("Invalid mutational context. All elements must have the same length, and the length must be 1, 3, or 5 characters.")
    }
  }
}


#Define the Transitions vs Transversions ratio --Ts_Tv_ratio {Ts:Tv} - for example 3:1 
if (any(grepl("--Ts_Tv_ratio", commandArgs()))) {
  ts_tv_ratio_pos <- which(grepl("--Ts_Tv_ratio", commandArgs()))
  
  if (length(ts_tv_ratio_pos) > 0 && ts_tv_ratio_pos < length(commandArgs())) {
    ts_tv_ratio <- commandArgs()[ts_tv_ratio_pos + 1]
    
    if (grepl(":", ts_tv_ratio)) {
      ts_tv_values <- as.numeric(unlist(strsplit(ts_tv_ratio, ":")))
      
      if (length(ts_tv_values) == 2 && !any(is.na(ts_tv_values))) {
        ts_value <- ts_tv_values[1]
        tv_value <- ts_tv_values[2]
        
        tot_ts_tv<-ts_value + tv_value
        ts_prob<-ts_value/tot_ts_tv
        tv_prob<-tv_value/tot_ts_tv
        variant_prob_table<-data.frame(
          Nucleotide=c("A","A","A","T","T","T","C","C","C","G","G","G"),
          Variant=c("T","C","G","A","C","G","G","A","T","C","A","T"),
          Probability=rep("")
        )
        variant_prob_table$Probability<-variant_prob_table$Probability<-ifelse(variant_prob_table$Nucleotide == "A" & variant_prob_table$Variant == "G", ts_prob,
                                                                               ifelse(variant_prob_table$Nucleotide == "A" & variant_prob_table$Variant == "C",tv_prob/2,
                                                                                      ifelse(variant_prob_table$Nucleotide == "A" & variant_prob_table$Variant == "T",tv_prob/2,
                                                                                             ifelse(variant_prob_table$Nucleotide == "T" & variant_prob_table$Variant == "A",tv_prob/2,
                                                                                                    ifelse(variant_prob_table$Nucleotide == "T" & variant_prob_table$Variant == "C",ts_prob,
                                                                                                           ifelse(variant_prob_table$Nucleotide == "T" & variant_prob_table$Variant == "G",tv_prob/2,
                                                                                                                  ifelse(variant_prob_table$Nucleotide == "C" & variant_prob_table$Variant == "G",tv_prob/2,
                                                                                                                         ifelse(variant_prob_table$Nucleotide == "C" & variant_prob_table$Variant == "A",tv_prob/2,
                                                                                                                                ifelse(variant_prob_table$Nucleotide == "C" & variant_prob_table$Variant == "T",ts_prob,
                                                                                                                                       ifelse(variant_prob_table$Nucleotide == "G" & variant_prob_table$Variant == "C",tv_prob/2,
                                                                                                                                              ifelse(variant_prob_table$Nucleotide == "G" & variant_prob_table$Variant == "A",ts_prob,tv_prob/2)))))))))))
        
        print(paste0("Transition vs Transversion ratio is ", ts_tv_ratio))
      } else {
        stop("Invalid format for --Ts_Tv_ratio. Expected format 'number:number'.")
      }
    } else {
      stop("Invalid format for --Ts_Tv_ratio. Expected format 'number:number'.")
    }
  }
} else if (any(grepl("^--variant_probability_table$", commandArgs()))) {
  #load probability table if it is given
  variant_prob_table_arg_index <- which(grepl("^--variant_probability_table$", commandArgs(trailingOnly = TRUE)))
  if (variant_prob_table_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    variant_prob_table_file <- commandArgs(trailingOnly = TRUE)[variant_prob_table_arg_index + 1]
    variant_prob_table_path <- file.path("/usr/src/app/pipeline/", variant_prob_table_file)
    if (file.exists(variant_prob_table_path)) {
      variant_prob_table <- read.csv(variant_prob_table_path)
      variant_prob_table[] <- lapply(variant_prob_table, function(x) if(is.character(x)) gsub(" ", "", x) else x)
      print("Table with the probability of a nucleotide mutating to another nucleotide is loaded")
    } else {
      stop("The specified varaint probability table not found")
    }
  } 
} else {
  print("table with variant mutation probability not provided, using default, each nucleotide can mutate to other 3 variants with equal probability")
  variant_prob_table<-data.frame(
    Nucleotide=c("A","A","A","T","T","T","C","C","C","G","G","G"),
    Variant=c("T","C","G","A","C","G","G","A","T","C","A","T"),
    Probability=rep(100/3)
  )
}

# Load the sbs table
if (any(grepl("^--sbs_signature_table$", commandArgs()))) {
  sbs_arg_index <- which(grepl("^--sbs_signature_table$", commandArgs(trailingOnly = TRUE)))
  
  if (sbs_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    sbs_file <- commandArgs(trailingOnly = TRUE)[sbs_arg_index + 1]
    sbs_path <- file.path(sbs_file)
    if (file.exists(sbs_path)) {
      sbs_table <- read.csv(sbs_path)
      print("SBS table loaded")
    } else {
      stop("The specified file does not exist.")
    }
  }
} else {
  print("The sbs_file is not specified")
}


##############################################################
#                                                            #
#         III. W H A T  I S  I N  T H E  S A M P L E ?       #
#                                                            #
##############################################################


#Define the number of copies in the sample
if (any(grepl("--nr_copies", commandArgs()))) {
  arg_pos <- which(grepl("--nr_copies", commandArgs()))
  
  if (length(arg_pos) > 0 && arg_pos < length(commandArgs())) {
    nr_copies <- as.integer(commandArgs()[arg_pos + 1])
    
    if (is.na(nr_copies) || nr_copies <= 0 || nr_copies != as.integer(nr_copies)) {
      stop("Invalid number of copies. It must be a positive integer.")
    }
    
    print(paste0("Number of copies set to ", nr_copies))
    nr_copies <- rep(nr_copies, times = nrow(fasta_content))
  } 
} else if(any(grepl("^--multifasta_copies$", commandArgs()))) {
  multifasta_copies_table_arg_index <- which(grepl("^--multifasta_copies$", commandArgs(trailingOnly = TRUE)))
  if (multifasta_copies_table_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    multifasta_copies_table_file <- commandArgs(trailingOnly = TRUE)[multifasta_copies_table_arg_index + 1]
    multifasta_copies_table_path <- file.path("../..", multifasta_copies_table_file)
    if (file.exists(multifasta_copies_table_path)) {
      multifasta_copies_table<-read.csv(multifasta_copies_table_path)
      nr_copies<-as.numeric(multifasta_copies_table$copy_number)
      print("Table with a copy number for each multifasta file is loaded")
    }else {
      stop("The specified table not found")
    } 
  } 
  
} else {
  print("Copy number not provided, using default - 1000 copies, if multifasta file provided, 1000 copies of each fasta sequence will be used")
  nr_copies <- 100
  nr_copies<-rep(nr_copies, time=nrow(fasta_content))
}


#################################################################################################
#                                                                                               #
#         IV.  D E F I N I N G  N E C E S S A R Y  F U N C T I O N S  &  V A R I A B L E S      #
#                                                                                               #
#################################################################################################
if (exists("mut_context")) {
  context_nchar<-nchar(mut_context)[1]
  if (context_nchar == 1) {
    SQL_column<-"nucleotide"
    SQL_column_rev<-"comp_nucleotide"
  } else if (context_nchar == 3) {
    SQL_column<-"tri_context"
    SQL_column_rev<-"tri_context_rev_comp"
  } else if (context_nchar == 5) {
    SQL_column<-"penta_context"
    SQL_column_rev<-"penta_context_rev_comp"
  } else {
    stop("Length of mutational context did not met the length requirement, 1, 3 or 5")
  }
} else {
  SQL_column<-"nucleotide"
  SQL_column_rev<-"comp_nucleotide"
}



#function to generate variants in mutations_table
generate_variants<-function(nucleotide) {
  probs<-filter(variant_prob_table,Nucleotide == nucleotide)
  return(sample(probs$Variant, 1, prob = probs$Probability))
}

# function ot get the total number of row
get_total_rows <- function(fasta_header) {
  conn <- dbConnect(SQLite(), dbname = paste0("/usr/src/app/pipeline/SQL_database/", fasta_header, ".sqlite"))
  total_rows_query <- sprintf("SELECT MAX(ROWID) as max_row_id FROM %s", fasta_header)
  total_rows_result <- dbGetQuery(conn, total_rows_query)
  dbDisconnect(conn)
  
  return(total_rows_result$max_row_id)
}


#Function for finding the rows that match the context
querying_SQL_with_context <- function(fasta_header, SQL_column, SQL_column_rev, mut_context, start_row, end_row) {
  sql_context <- paste0("'", mut_context, "'", collapse = ", ")
  
  conn <- dbConnect(SQLite(), dbname = paste0("/usr/src/app/pipeline/SQL_database/", fasta_header, ".sqlite"))
  
  query_5 <- sprintf("SELECT ROWID as row_number FROM %s WHERE %s IN (%s) AND ROWID BETWEEN %d AND %d", 
                     fasta_header, SQL_column, sql_context, start_row, end_row)
  result_5 <- dbGetQuery(conn, query_5)
  
  mut_positions_5 <- result_5$row_number
  
  query_3 <- sprintf("SELECT ROWID as row_number FROM %s WHERE %s IN (%s) AND ROWID BETWEEN %d AND %d", 
                     fasta_header, SQL_column_rev, sql_context, start_row, end_row)
  result_3 <- dbGetQuery(conn, query_3)
  mut_positions_3 <- result_3$row_number
  
  dbDisconnect(conn)
  
  # Return the positions in a list
  return(list(mut_positions_5 = mut_positions_5, mut_positions_3 = mut_positions_3))
}

#Function that subselects the SQL database table, adds column variant and mutates it to desired variant

get_mutations_table <- function(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants, start_row, end_row) {
  conn <- dbConnect(SQLite(), dbname = paste0("/usr/src/app/pipeline/SQL_database/", fasta_header, ".sqlite"))
  
  selected_positions_5 <- mut_positions_5_selected[mut_positions_5_selected >= start_row & mut_positions_5_selected <= end_row]
  selected_positions_5<-sort(selected_positions_5)
  selected_positions_5<-unique(selected_positions_5)
  selected_positions_3 <- mut_positions_3_selected[mut_positions_3_selected >= start_row & mut_positions_3_selected <= end_row]
  selected_positions_3<-sort(selected_positions_3)
  selected_positions_3<-unique(selected_positions_3)

  mutations_table_5 <- data.frame()
  mutations_table_3 <- data.frame()
  
  # Query the 5' strand
  if (length(selected_positions_5) > 0) {
    row_ids_str_5 <- paste(selected_positions_5, collapse = ", ")
    query_5 <- sprintf("SELECT * FROM %s WHERE ROWID IN (%s)", fasta_header, row_ids_str_5)
    mutations_table_5 <- dbGetQuery(conn, query_5)
    mutations_table_5$position <- selected_positions_5
    mutations_table_5$Variant <- sapply(mutations_table_5$nucleotide, function(nuc) {
      ifelse(nuc == "N", "N", generate_variants(nuc))
    })
  }
  
  # Query the 3' strand
  if (length(selected_positions_3) > 0) {
    row_ids_str_3 <- paste(selected_positions_3, collapse = ", ")
    query_3 <- sprintf("SELECT * FROM %s WHERE ROWID IN (%s)", fasta_header, row_ids_str_3)
    mutations_table_3 <- dbGetQuery(conn, query_3)
    mutations_table_3$position <- selected_positions_3
    mutations_table_3$Variant <- sapply(mutations_table_3$nucleotide, function(nuc) {
      ifelse(nuc == "N", "N", generate_variants(nuc))
    })
    mutations_table_3$Variant <- sapply(mutations_table_3$Variant, function(var) {
      ifelse(var == "N", "N", chartr("ACGT", "TGCA", var))
    })
    
  }
  
  dbDisconnect(conn)
  
  return(rbind(mutations_table_5, mutations_table_3))
}




#######################################################################################################
#                                                                                                     #
#         V.  M U T A T E D  G E N O M E   F R A C T I O N   O R   M U T A T I O N   R A T E ?        #
#                                                                                                     #
#######################################################################################################

setwd(output_folder) 

if (any(grepl("--mut_genome_fraction", commandArgs()))) {
  arg_pos <- which(grepl("--mut_genome_fraction", commandArgs()))
  
  if (length(arg_pos) > 0 && arg_pos < length(commandArgs())) {
    mut_genome <- as.numeric(commandArgs()[arg_pos + 1])
    if (mut_genome <= 0 || mut_genome > 1) {
      stop("mut_genome must be greater than 0 and less or equal to 1.")
    }
    print(paste0("Fraction of all genomes that will get a mutation is set to ", mut_genome))
  }
  
  # Check for fraction of positions
  if (any(grepl("--fraction_positions", commandArgs()))) {
    arg_pos <- which(grepl("--fraction_positions", commandArgs()))
    if (length(arg_pos) > 0 && arg_pos < length(commandArgs())) {
      pos_fraction <- as.numeric(commandArgs()[arg_pos + 1])
      if (pos_fraction <= 0 || pos_fraction > 1) {
        stop("pos_fraction must be greater than 0 and less or equal to 1.")
      }
      print(paste0("Fraction of positions that will be mutated is ", pos_fraction))
    } 
  } else {
    print("Fraction of mutation not provided, using the default of 0.05 (5%)")
    pos_fraction <- 0.05
  }
  
} else if (any(grepl("--specific_mut_rate", commandArgs()))) {
  print("Fraction of mutated genome not provided, using specified mutation rate")
  arg_pos <- which(grepl("--specific_mut_rate", commandArgs()))
  if (length(arg_pos) > 0 && arg_pos < length(commandArgs())) {
    specific_mut_rate <- as.numeric(commandArgs()[arg_pos + 1])
    if (is.na(specific_mut_rate) || specific_mut_rate < 0 || specific_mut_rate > 1) {
      stop("Invalid specific mutation rate. It must be a positive number.")
    }
    print(paste0("Specific mutation rate:", specific_mut_rate))
  } else {
    print("Specific mutation rate not defined, continuing with mutated genome fraction.")
    specific_mut_rate <- NULL
  }
  
} else {
  stop("None of the genome mutation options was provided")
}

################################################################################
#                                                                              #
#           VA    M U T A T E D   G E N O M E   F R A C T I O N                #     
#                                                                              #
################################################################################


print("Multiplying and mutating genomes in provided fasta file")


lapply(1:nrow(fasta_content), function(seq_index) {
  
  fasta_header <- fasta_content$fasta_names[seq_index]
  
  bed_data_filtered<-bed_data[bed_data$V1 == fasta_header,]
  
  # split filtered bed file into processing batches
  
  split_into_batches <- function(data, num_cores) {
    split(data, cut(seq(nrow(data)), num_cores, labels = FALSE))
  }
  
  batches <- split_into_batches(bed_data_filtered, num_cores)
  
  
  #if mut_genome (coming from mut_genome_fraction command) exisits and sbs_table doesn't exist
  if (exists("mut_genome") && !exists("sbs_table")) {
    # and if mut_context is given
    if (exists("mut_context")) {

      process_batch <- function(batch) {
        apply(batch, 1, function(row) {
          start_pos = as.numeric(row[2])
          end_pos = as.numeric(row[3])
          
          querying_SQL_with_context(fasta_header, SQL_column, SQL_column_rev, mut_context, start_pos, end_pos)
        })
      }
      
      results <- mclapply(batches, process_batch, mc.cores = num_cores)
      
      mut_positions_5 <- c()
      mut_positions_3 <- c()
      
      extract_positions <- function(results) {
        # Loop through each batch
        for (batch_key in names(results)) {
          batch = results[[batch_key]]
          
          for (sub_batch_key in names(batch)) {
            sub_batch = batch[[sub_batch_key]]
            
            if (!is.null(sub_batch$mut_positions_5)) {
              mut_positions_5 <<- c(mut_positions_5, sub_batch$mut_positions_5)
            }
            
            if (!is.null(sub_batch$mut_positions_3)) {
              mut_positions_3 <<- c(mut_positions_3, sub_batch$mut_positions_3)
            }
          }
        }
      }
      
      extract_positions(results)
      
      mut_positions <- sort(unique(c(mut_positions_5, mut_positions_3)))
      
    } else {
      print("Trinucleotide context is not provided, using all positions")
      
      mut_positions <- unlist(lapply(1:nrow(bed_data_filtered), function(i) seq(bed_data_filtered$V2[i], bed_data_filtered$V3[i])))
    }
    
    
    tot_mut_number<-round(length(mut_positions)*pos_fraction)
    
    print(paste0("Randomly picking ", tot_mut_number," positions to be mutated"))
    
    picked_mut_positions <- sample(mut_positions, tot_mut_number,replace = FALSE)
    picked_mut_positions <- sort(picked_mut_positions)
    
    print("Making a table that stores the position, orignial nucleotide, tri-nucleotide context, and to which nucleotide it will be mutated")
    
    if (exists("mut_positions_5")) {
      
      
      mut_positions_5_selected<-picked_mut_positions[picked_mut_positions %in% mut_positions_5]
      mut_positions_3_selected<-picked_mut_positions[picked_mut_positions %in% mut_positions_3]
      
      worker_seeds <- NULL
      
      if (exists("seed_value")) {
        # Generate worker seeds since the main seed is set
        worker_seeds <- sample(1:1e9, num_cores, replace = FALSE)
      }
      
      process_mutation_batch <- function(batch, seed) {
        set.seed(seed)
        results <- apply(batch, 1, function(row) {
          start_pos <- as.numeric(row[2]) 
          end_pos <- as.numeric(row[3]) 
          
          get_mutations_table(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants, start_pos, end_pos)
        })
        
        do.call(rbind, results)
      }
      
      mutations_table_list <- mclapply(seq_along(batches), function(i) {
        process_mutation_batch(batches[[i]], worker_seeds[i])
      }, mc.cores = num_cores)
      
      mutations_table <- do.call(rbind, mutations_table_list)
      
      
      
    } else { 
      #if the mut_context is not provided
      mut_positions_5_selected<-sample(picked_mut_positions, size = round(0.5 * length(picked_mut_positions)), replace = FALSE)
      mut_positions_3_selected <- setdiff(picked_mut_positions, mut_positions_5_selected)
      
      
      worker_seeds <- NULL
      
      if (exists("seed_value")) {
        worker_seeds <- sample(1:1e9, num_cores, replace = FALSE)
      }
      
      process_mutation_batch <- function(batch, seed) {
        set.seed(seed)
        # Process each row in the batch data frame
        results <- apply(batch, 1, function(row) {
          start_pos <- as.numeric(row[2]) 
          end_pos <- as.numeric(row[3]) 
          
          get_mutations_table(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants, start_pos, end_pos)
        })
        
        do.call(rbind, results)
      }
      
      mutations_table_list <- mclapply(seq_along(batches), function(i) {
        process_mutation_batch(batches[[i]], worker_seeds[i])
      }, mc.cores = num_cores)
      
      mutations_table <- do.call(rbind, mutations_table_list)
      
    }
    # Write mutations table
    mutations_table<-mutations_table[order(mutations_table$position),]
    mutations_table_name<-paste0(fasta_header,"_mutations_table.csv")
    fwrite(mutations_table, mutations_table_name, sep = ";")
    
  
    
    MGs_number<-round(nr_copies[seq_index]*mut_genome)
    
    sample_table_NMG<-data.frame(
      genome_name=paste0(fasta_header, "_NMG"),
      variant_position = "No_mutations",
      copy_number = nr_copies[seq_index] - MGs_number 
    )
    NMG_table_temp_name<-paste0(fasta_header, "_NMG_temp.csv")
    fwrite(sample_table_NMG, NMG_table_temp_name, sep = ";")
    

    worker_seeds <- NULL
    
    if (exists("seed_value")) {
      worker_seeds <- sample(1:1e9, MGs_number, replace = FALSE)
    } 
    
    
    mutations_table$pos_var<-paste0(mutations_table$position, mutations_table$Variant)
    all_possible_mutations<-mutations_table$pos_var
    
    output_files <- mclapply(1:MGs_number, function(i) {
      set.seed(worker_seeds[i])
      num_mutations <- sample(length(all_possible_mutations), 1)
      selected_mutations <- sample(all_possible_mutations, num_mutations, replace = FALSE)
      selected_mutations <- selected_mutations[order(as.numeric(gsub("[A-Za-z]", "", selected_mutations)))]
      
      selected_mutations <- paste(selected_mutations, collapse = ",")
      
      temp_file <- paste0("temp_rows_MG_", Sys.getpid(), "_", i, ".csv")
      fwrite(data.table(genome_name = fasta_header, variant_position = selected_mutations, copy_number = 1),
             temp_file, sep = ";", quote = FALSE, row.names = FALSE, col.names = TRUE)
      
      temp_file
    }, mc.cores = num_cores)
    
    
    MG_table_temp_name <- paste0(fasta_header,"_MG_temp.csv")
    if(file.exists(MG_table_temp_name)) file.remove(MG_table_temp_name)
    
    system(paste("head -n 1", output_files[1], ">", MG_table_temp_name))
    
    if(length(output_files) > 1) {
      system(paste("tail -n +2 -q", paste(output_files, collapse = " "), ">>", MG_table_temp_name))
    }
    
    output_files<-paste0("./",unlist(output_files))
    
    file.remove(output_files)
    
    awk_cmd <- paste(
      "awk -F';' 'BEGIN {OFS=\";\"} {if(NR>1) $1=$1\"_MG\"(NR-1); print}'",
      MG_table_temp_name, "> temp_file && mv temp_file", MG_table_temp_name
    )
    system(awk_cmd)
    

    sample_table_name<-paste0(fasta_header, "_sample_table.csv")
    
    header_cmd <- paste("head -n 1", NMG_table_temp_name, ">" , sample_table_name)
    system(header_cmd)
    
    tail_cmd_1 <- paste("tail -n +2", NMG_table_temp_name, ">>", sample_table_name)
    system(tail_cmd_1)
    
    tail_cmd_2 <- paste("tail -n +2", MG_table_temp_name, ">>" , sample_table_name)
    system(tail_cmd_2)
    
    temp_NMG_MG<-paste0("./", c(NMG_table_temp_name, MG_table_temp_name))
    file.remove(temp_NMG_MG)

    
  } else if (exists("specific_mut_rate")) {
    
    ################################################################################
    #                                                                              #
    #           VB   S P E C I F I C   M U T A T I O N   R A T E                   #     
    #                                                                              #
    ################################################################################
    

    print("mutating with specific mut rate")
    
    interval_lengths <- bed_data_filtered$V3 - bed_data_filtered$V2 + 1
    
    total_exome_length <- sum(interval_lengths)
    
    lambda<-specific_mut_rate * total_exome_length
    
    prob_MGs<-rpois(n=nr_copies[seq_index], lambda = lambda)
    
    has_mutation<-prob_MGs > 0
    n_mutations<-prob_MGs[has_mutation]
    
    if(sum(n_mutations) > 0){
      
      MGs_number<-length(n_mutations)
      
      sample_table_NMG<-data.frame(
        genome_name=paste0(fasta_header, "_NMG"),
        variant_position = "No_mutations",
        copy_number = nr_copies[seq_index] - MGs_number 
      )

      MGs_table<-data.frame(
        genome_name = paste0(fasta_header, "_MG", 1:MGs_number),
        nr_mutations = n_mutations,
        positions = rep("", times = MGs_number),
        variant_position = rep("", times = MGs_number)
      )
      
      if (exists("mut_context")) {

        process_batch <- function(batch) {
          apply(batch, 1, function(row) {
            # Convert start and end positions to numeric
            start_pos = as.numeric(row[2])
            end_pos = as.numeric(row[3])
            
            querying_SQL_with_context(fasta_header, SQL_column, SQL_column_rev, mut_context, start_pos, end_pos)
          })
        }
        
        results <- mclapply(batches, process_batch,mc.cores = num_cores)
        
        mut_positions_5 <- c()
        mut_positions_3 <- c()
        
        extract_positions <- function(results) {
          for (batch_key in names(results)) {
            batch = results[[batch_key]]
            
            for (sub_batch_key in names(batch)) {
              sub_batch = batch[[sub_batch_key]]
              
              if (!is.null(sub_batch$mut_positions_5)) {
                mut_positions_5 <<- c(mut_positions_5, sub_batch$mut_positions_5)
              }
              
              if (!is.null(sub_batch$mut_positions_3)) {
                mut_positions_3 <<- c(mut_positions_3, sub_batch$mut_positions_3)
              }
            }
          }
        }
        
        extract_positions(results)
        
        mut_positions <- sort(unique(c(mut_positions_5, mut_positions_3))) 
        
        MGs_table$positions <- sapply(MGs_table$nr_mutations, function(x) {
          paste(sample(mut_positions, x, replace = F), collapse = ",")
        })
        
        picked_mut_positions<-unique(sort(as.numeric(unlist(strsplit(MGs_table$positions, ",")))), decreasing = FALSE)

        mut_positions_5_selected<-picked_mut_positions[picked_mut_positions %in% mut_positions_5]
        mut_positions_3_selected<-picked_mut_positions[picked_mut_positions %in% mut_positions_3]
        
        worker_seeds <- NULL
        
        if (exists("seed_value")) {
          worker_seeds <- sample(1:1e9, num_cores, replace = FALSE)
        }
        
        process_mutation_batch <- function(batch, seed) {
          set.seed(seed)
          results <- apply(batch, 1, function(row) {
            start_pos <- as.numeric(row[2])  
            end_pos <- as.numeric(row[3]) 
            
            get_mutations_table(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants, start_pos, end_pos)
          })
          
          do.call(rbind, results)
        }
        
        mutations_table_list <- mclapply(seq_along(batches), function(i) {
          process_mutation_batch(batches[[i]], worker_seeds[i])
        }, mc.cores = num_cores)
        
        mutations_table <- do.call(rbind, mutations_table_list)
        
        
      } else {
        #when mutational context is not given every possible position can mutate
        mut_positions <- unlist(lapply(1:nrow(bed_data_filtered), function(i) seq(bed_data_filtered$V2[i], bed_data_filtered$V3[i])))
        MGs_table$positions <- sapply(MGs_table$nr_mutations, function(x) {
          paste(sample(mut_positions, x, replace = F), collapse = ",")
        })
        picked_mut_positions<-paste(unique(sort(as.numeric(unlist(strsplit(MGs_table$positions, ",")))), decreasing = FALSE), collapse = ",")
        
        mut_positions_5_selected<-sample(picked_mut_positions, size = round(0.5 * length(picked_mut_positions)), replace = FALSE)
        mut_positions_3_selected <- as.numeric(unlist(strsplit(setdiff(picked_mut_positions, mut_positions_5_selected), ",")))
        
        worker_seeds <- NULL
        
        if (exists("seed_value")) {
          worker_seeds <- sample(1:1e9, num_cores, replace = FALSE)
        }
        
        process_mutation_batch <- function(batch, seed) {
          set.seed(seed)
          results <- apply(batch, 1, function(row) {
            start_pos <- as.numeric(row[2])  
            end_pos <- as.numeric(row[3]) 
            
            get_mutations_table(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants, start_pos, end_pos)
          })
          
          do.call(rbind, results)
        }
        
        mutations_table_list <- mclapply(seq_along(batches), function(i) {
          process_mutation_batch(batches[[i]], worker_seeds[i])
        }, mc.cores = num_cores)
        
        mutations_table <- do.call(rbind, mutations_table_list)
      }
      
      setDT(MGs_table)
      setDT(mutations_table)
      
      mutations_table[, position := as.character(position)]
      get_variants <- function(pos_str) {
        pos_vec <- unlist(strsplit(pos_str, ","))
        variants <- mutations_table[position %in% pos_vec, .(position, Variant)][order(position)]
        paste0(variants$position, variants$Variant, collapse = ",")
      }
      
      MGs_table[, variant_position := sapply(positions, get_variants)]
      MGs_table<-as.data.frame(MGs_table)
      
      
      sample_table_MG<-MGs_table %>% 
        select(genome_name, variant_position) %>% 
        mutate(copy_number = 1) %>%
        group_by(variant_position) %>%
        summarise(genome_name = dplyr::first(genome_name), copy_number = n()) %>%
        ungroup() %>%  
        select(genome_name, variant_position, copy_number) 
      
      
      sample_table<-rbind(sample_table_NMG, sample_table_MG)
      sample_table<-sample_table[sample_table$copy_number !=0,]
      
      sample_table_name<-paste0(fasta_header, "_sample_table.csv")
      fwrite(sample_table, sample_table_name, sep = ";")
      
      mutations_table_name<-paste0(fasta_header,"_mutations_table.csv")
      fwrite(mutations_table, mutations_table_name, sep = ";")
      
      print("*_mutations_table.csv contains the information about mutations that was insreted into mutated genoms (one file per sequence in fasta file)")
      print("*_sample_table.csv contains the information about how many genomes got a mutation and which positions/nucleotide were mutated")
      
      
    } else {
      print(paste("No copies of genome",fasta_header,"were mutated. The sample table with only non-mutated genomes will be made for this genome, copy number will correspond to the number of copies"))
      sample_table<-data.frame(
        genome_name=paste0(fasta_header, "_NMG"),
        variant_position = "No_mutations",
        # copy number is equal to total number of copies minus the fraction of sample that will be mutated
        copy_number = nr_copies[seq_index]
      )
      
      sample_table_name<-paste0(fasta_header, "_sample_table.csv")
      fwrite(sample_table, sample_table_name, sep = ";")
      
    }
    
    
    
  } else if (exists("sbs_table") && exists("mut_genome")) {
    
    ################################################################################################################
    #                                                                                                              #
    #           VC R E P L I C A T I N G  S I G N A T U R E S  F R O M  S B S  T A B L E S                         #     
    #                                                                                                              #
    ################################################################################################################
    
    print("Mutatiting genomes using sbs signtaure table")

    transformed_data <- sbs_table %>%
      mutate(
        context = str_extract(X, "^[A-Z]") %>% paste0(str_extract(X, "(?<=\\[)[A-Z](?=>)"), str_extract(X, "(?<=\\])[A-Z]")),
        variant = str_extract(X, "(?<=\\>)[A-Z]")
      ) %>%
      select(3,4,2) 
    
    names(transformed_data)[3] <- "occurance_prob"

    contexts_sum<-transformed_data %>% group_by(context) %>% summarise(occurance_prob = sum(occurance_prob))
    
    variant_probability <- transformed_data %>%
      pivot_wider(names_from = variant, values_from = occurance_prob)
    
    variant_probability[is.na(variant_probability)] <- 0
    variant_probability<-variant_probability%>%
      rowwise() %>%
      mutate(
        sum = sum(c_across(A:C)),
        A = A / sum,
        G = G / sum,
        T = T / sum,
        C = C / sum
      ) %>%
      select(-sum)
    
    conn_str <- paste0("/usr/src/app/pipeline/SQL_database/", fasta_header, ".sqlite")
    
    process_batch <- function(batch_data, conn_str, contexts, fasta_header) {
      con <- dbConnect(RSQLite::SQLite(), conn_str)
      on.exit(dbDisconnect(con))
      
      temp_results <- data.frame(context = contexts, genome_occurence = rep(0, length(contexts)))
      
      for (i in 1:nrow(batch_data)) {
        exon_start <- batch_data[i, 2]
        exon_end <- batch_data[i, 3]
        
        query <- sprintf(
          "SELECT tri_context, tri_context_rev_comp, COUNT(*) as count FROM %s WHERE (tri_context IN (%s) OR tri_context_rev_comp IN (%s)) AND ROWID BETWEEN %d AND %d GROUP BY tri_context, tri_context_rev_comp",
          fasta_header, paste0("'", contexts, "'", collapse = ","), paste0("'", contexts, "'", collapse = ","), exon_start, exon_end
        )
        
        results <- dbGetQuery(con, query)
        
        for (j in 1:nrow(results)) {
          context <- results$tri_context[j]
          count <- results$count[j]
          temp_results$genome_occurence[temp_results$context == context] <- temp_results$genome_occurence[temp_results$context == context] + count
          
          context_rev <- results$tri_context_rev_comp[j]
          count_rev <- results$count[j]
          temp_results$genome_occurence[temp_results$context == context_rev] <- temp_results$genome_occurence[temp_results$context == context_rev] + count_rev
        }
      }
      
      return(temp_results)
    }
    context_genome_occurence<-contexts_sum %>% select(-occurance_prob) %>% mutate(genome_occurence = 0)
    contexts <- context_genome_occurence$context
    
    results_list <- mclapply(batches, function(batch_data) {
      process_batch(batch_data, conn_str, contexts, fasta_header)
    }, mc.cores = num_cores)
    
    combined_results <- Reduce(function(x, y) {
      x$genome_occurence <- x$genome_occurence + y$genome_occurence
      return(x)
    }, results_list)
    
    context_genome_occurence$genome_occurence <- combined_results$genome_occurence
    

    context_genome_occurence <- context_genome_occurence %>% arrange(context)
    contexts_sum <- contexts_sum %>% arrange(context)
    
    context_genome_occurence$prob<-context_genome_occurence$genome_occurence/sum(context_genome_occurence$genome_occurence)
    
    contexts_sum$adjusted<-contexts_sum$occurance_prob*context_genome_occurence$prob
    
    contexts_sum$adjusted<-contexts_sum$occurance_prob*context_genome_occurence$prob
    
    total_occurrences <- sum(context_genome_occurence$genome_occurence)
    

    tot_mut_number<-round(total_occurrences*pos_fraction)
    
    probabilities <- contexts_sum$adjusted
    
    mutations <- rmultinom(1, tot_mut_number, probabilities)
    
    mutations_df <- data.frame(
      context = contexts_sum$context,
      mutations = as.vector(mutations)
    )
    
    if (num_cores > 32) {
      num_cores_context <- 32
    } else {num_cores_context<-num_cores}
    
    worker_seeds <- NULL
    
    if (exists("seed_value")) {
      worker_seeds <- sample(1:1e9, num_cores, replace = FALSE)
    }
    
    conn <- dbConnect(SQLite(), dbname = conn_str)
    
    process_context <- function(i, mutations_df, bed_data_filtered, conn_str, fasta_header, worker_seeds) {
      set.seed(worker_seeds[i %% length(worker_seeds) + 1])
      
      context <- mutations_df$context[i]
      mutations <- mutations_df$mutations[i]
      
      all_positions_context <- c()
      batch_size <- 100  
      
      conn <- dbConnect(SQLite(), dbname = conn_str)
      on.exit(dbDisconnect(conn), add = TRUE) 
      
      bed_batches <- split(bed_data_filtered, ceiling(seq_along(bed_data_filtered$V2) / batch_size))
      
      for (bed_batch in bed_batches) {
        where_clauses <- apply(bed_batch, 1, function(row) {
          sprintf("(ROWID BETWEEN %s AND %s)", row["V2"], row["V3"])
        })
        combined_where_clause <- paste(where_clauses, collapse = " OR ")
        
        query <- sprintf(
          "SELECT ROWID FROM %s WHERE (tri_context = '%s' OR tri_context_rev_comp = '%s') AND (%s)",
          fasta_header, context, context, combined_where_clause
        )
        
        result <- dbGetQuery(conn, query)$rowid
        all_positions_context <- c(all_positions_context, result)
      }
      
      if (length(all_positions_context) > 0) {
        sampled_positions <- sample(all_positions_context, size = mutations)
        return(list(context = context, sampled_positions = sampled_positions))
      } else {
        return(NULL)
      }
    }
    
    results <- mclapply(1:nrow(mutations_df), function(i) {
      process_context(i, mutations_df, bed_data_filtered, conn_str, fasta_header, worker_seeds)
    }, mc.cores = num_cores_context)
    
    all_sampled_positions <- unlist(lapply(results, function(x) x$sampled_positions))
    
    dbDisconnect(conn)
    
    get_mutations_table_sbs <- function(fasta_header, mut_positions_selected, start_row, end_row) {
      conn <- dbConnect(SQLite(), dbname = conn_str)
      
      selected_positions <- mut_positions_selected[mut_positions_selected >= start_row & mut_positions_selected <= end_row]
      selected_positions<-sort(selected_positions)
      
      mutations_table <- data.frame()
      
      
      if (length(selected_positions) > 0) {
        row_ids <- paste(selected_positions, collapse = ", ")
        query <- sprintf("SELECT * FROM %s WHERE ROWID IN (%s)", fasta_header, row_ids)
        mutations_table <- dbGetQuery(conn, query)
        mutations_table$position <- selected_positions
        
      }
      
      dbDisconnect(conn)
      
      return(rbind(mutations_table))
    }
    
    
    process_mutation_batch <- function(batch) {
      results <- apply(batch, 1, function(row) {
        start_pos <- as.numeric(row[2])  
        end_pos <- as.numeric(row[3])    
        
        get_mutations_table_sbs(fasta_header, all_sampled_positions, start_pos, end_pos)
      })
      
      do.call(rbind, results)
    }
    
    mutations_table_list <- mclapply(batches, process_mutation_batch, mc.cores = num_cores)
    
    mutations_table <- do.call(rbind, mutations_table_list)
    
    generate_variant <- function(probabilities) {
      bases <- c("A", "G", "T", "C")
      sample(bases, 1, prob = probabilities)
    }
    
    get_complement <- function(nucleotide) {
      complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
      return(complement[nucleotide])
    }
    
    process_row <- function(index) {
      if (!is.null(worker_seeds)) {
        set.seed(worker_seeds[(index %% num_cores) + 1])
      }
      
      row <- mutations_table[index, ]
      tri_context <- row$tri_context
      tri_context_rev_comp <- row$tri_context_rev_comp
      
      if (tri_context %in% variant_probability$context) {
        probabilities <- variant_probability %>% filter(context == tri_context) %>% select(A, G, T, C)
        generated_variant <- generate_variant(as.numeric(probabilities))
      } else if (tri_context_rev_comp %in% variant_probability$context) {
        probabilities <- variant_probability %>% filter(context == tri_context_rev_comp) %>% select(A, G, T, C)
        generated_variant <- generate_variant(as.numeric(probabilities))
        generated_variant <- get_complement(generated_variant)
      } else {
        generated_variant <- NA
      }
      
      row$Variant <- generated_variant
      return(row)
    }
    
    results <- mclapply(1:nrow(mutations_table), function(i) {
      process_row(i)
    }, mc.cores = num_cores)
    
    results_df <- do.call(rbind, results)
    
    mutations_table <- as_tibble(results_df)
    
    
    mutations_table<-mutations_table[order(mutations_table$position),]
    mutations_table_name<-paste0(fasta_header,"_mutations_table.csv")
    fwrite(mutations_table, mutations_table_name, sep = ";")
    
    MGs_number<-round(nr_copies[seq_index]*mut_genome)
    
    sample_table_NMG<-data.frame(
      genome_name=paste0(fasta_header, "_NMG"),
      variant_position = "No_mutations",
      copy_number = nr_copies[seq_index] - MGs_number 
    )
    NMG_table_temp_name<-paste0(fasta_header, "_NMG_temp.csv")
    fwrite(sample_table_NMG, NMG_table_temp_name, sep = ";")
    

    worker_seeds <- NULL
    
    if (exists("seed_value")) {
      worker_seeds <- sample(1:1e9, MGs_number, replace = FALSE)
    } 
    
    
    mutations_table$pos_var<-paste0(mutations_table$position, mutations_table$Variant)
    all_possible_mutations<-mutations_table$pos_var
    
    output_files <- mclapply(1:MGs_number, function(i) {
      set.seed(worker_seeds[i])
      num_mutations <- sample(length(all_possible_mutations), 1)
      selected_mutations <- sample(all_possible_mutations, num_mutations, replace = FALSE)
      selected_mutations <- selected_mutations[order(as.numeric(gsub("[A-Za-z]", "", selected_mutations)))]
      
      selected_mutations <- paste(selected_mutations, collapse = ",")
      
      temp_file <- paste0("temp_rows_MG_", Sys.getpid(), "_", i, ".csv")
      fwrite(data.table(genome_name = fasta_header, variant_position = selected_mutations, copy_number = 1),
             temp_file, sep = ";", quote = FALSE, row.names = FALSE, col.names = TRUE)
      
      temp_file
    }, mc.cores = num_cores)
    
    
    MG_table_temp_name <- paste0(fasta_header,"_MG_temp.csv")
    if(file.exists(MG_table_temp_name)) file.remove(MG_table_temp_name)
    
    system(paste("head -n 1", output_files[1], ">", MG_table_temp_name))
    
    if(length(output_files) > 1) {
      system(paste("tail -n +2 -q", paste(output_files, collapse = " "), ">>", MG_table_temp_name))
    }
    
    output_files<-paste0("./",unlist(output_files))
    
    file.remove(output_files)
    
    awk_cmd <- paste(
      "awk -F';' 'BEGIN {OFS=\";\"} {if(NR>1) $1=$1\"_MG\"(NR-1); print}'",
      MG_table_temp_name, "> temp_file && mv temp_file", MG_table_temp_name
    )
    system(awk_cmd)
    
    sample_table_name<-paste0(fasta_header, "_sample_table.csv")
    
    header_cmd <- paste("head -n 1", NMG_table_temp_name, ">" , sample_table_name)
    system(header_cmd)
    
    tail_cmd_1 <- paste("tail -n +2", NMG_table_temp_name, ">>", sample_table_name)
    system(tail_cmd_1)
    
    tail_cmd_2 <- paste("tail -n +2", MG_table_temp_name, ">>" , sample_table_name)
    system(tail_cmd_2)
    
    temp_NMG_MG<-paste0("./", c(NMG_table_temp_name, MG_table_temp_name))
    file.remove(temp_NMG_MG)
    
    
  } })


print("Process finished.")


