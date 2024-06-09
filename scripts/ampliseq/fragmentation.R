library(Biostrings)
library(foreach)
library(iterators)
library(doParallel)
library(data.table)
library(Rcpp)



##################################################
#                                                #
#   I.   D E F I N I N G   A R G U M E N T S     #
#                                                #
##################################################


#Setting the working direcotry
if (any(grepl("--input_directory", commandArgs()))) { 
  input_arg <- which(grepl("--input_directory", commandArgs())) 
  if (length(input_arg) > 0 && input_arg < length(commandArgs())) { 
    input_folder <- commandArgs()[input_arg + 1] 
  } 
} 

setwd(input_folder)

if (any(grepl("--output_directory", commandArgs()))) { 
  output_arg <- which(grepl("--output_directory", commandArgs())) 
  if (length(output_arg) > 0 && output_arg < length(commandArgs())) { 
    output_folder <- commandArgs()[output_arg + 1] 
  } 
}


# argument --set_seed {integer} - sets a seed for reproducibility of random processes 
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

##Define fraction of a sample that will go to fragmentation
if (any(grepl("--fragmentation_fraction", commandArgs()))) {
  fragmentation_fraction_pos <- which(grepl("--fragmentation_fraction", commandArgs()))
  
  if (length(fragmentation_fraction_pos) > 0 && fragmentation_fraction_pos < length(commandArgs())) {
    fragmentation_fraction <- as.numeric(commandArgs()[fragmentation_fraction_pos + 1])
    
    # Check if fragmentation_fraction is between 0 and 1
    if (is.na(fragmentation_fraction) || fragmentation_fraction < 0 || fragmentation_fraction > 1) {
      stop("Invalid fragmentation fraction. The value must be between 0 and 1, hope that you learned you lesson...")
    }
    
    print(paste0("Fragmentation fraction is set to: ", fragmentation_fraction, "hope that you know what you are doing..."))
  }
} else {
  print("Fragmentation fraction not defined, we will use default 1 (100%)")
  fragmentation_fraction <-1
}

## Is genome circular or not
if (any(grepl("--circular", commandArgs()))) {
  circular_pos <- which(grepl("--circular", commandArgs()))
  if (length(circular_pos) > 0 && circular_pos < length(commandArgs())) {
    circular_option <- commandArgs()[circular_pos + 1]
    if (circular_option == "Y") {
      print("Genome specifed as circular")
    } else { 
      print("Genome specifed as linear")}
  } 
} else{
  circular_option<- "Y"
  print("Genome was neither specifed as circular nor linear, using defaul - circular")
}

##Define the length of a fragment that will be produced during fragmentation
if (any(grepl("--fragment_length", commandArgs()))) {
  fragment_length_pos <- which(grepl("--fragment_length", commandArgs()))
  
  if (length(fragment_length_pos) > 0 && fragment_length_pos < length(commandArgs())) {
    fragment_length_arg <- commandArgs()[fragment_length_pos + 1]
    
    if (grepl(":", fragment_length_arg)) {
      print("Fragment length defined as a range")
      
      # Split the argument and convert to numeric
      range_values <- as.integer(unlist(strsplit(fragment_length_arg, ":")))
      
      # Check if range_values are integers and greater than zero
      if (length(range_values) == 2 && all(!is.na(range_values)) && all(range_values > 0)) {
        # Check if the range is valid (start <= end)
        if (range_values[1] <= range_values[2]) {
          fragment_length <- range_values[1]:range_values[2]
          print(paste0("Fragment length range: ", min(fragment_length),"-",max(fragment_length)))
        } else {
          stop("Invalid range")
        }
      } else {
        stop("Invalid range.")
      }
    } else {
      # If the argument is a single value, convert it to numeric
      fragment_length <- as.numeric(fragment_length_arg)
      print(paste("Fragment length defined as a single value", fragment_length))
      if (is.na(fragment_length) || fragment_length < 0) {
        stop("Invalid fragment length")
      }
    }
  }
} else {
  print("Fragment length values not provided, using deafult range 250:450")
  fragment_length<-250:450
}

# argument --cores {integer} defines number of cores, if not given process is run on 1 core sequentially
if (any(grepl("--cores", commandArgs()))) {
  core_arg <- which(grepl("--cores", commandArgs()))
  if (length(core_arg) > 0 && core_arg < length(commandArgs())) {
    num_cores <- as.integer(commandArgs()[core_arg + 1])
    
    # Check if num_cores is a positive integer
    if (is.na(num_cores) || num_cores <= 0 || num_cores != floor(num_cores)) {
      stop("Invalid number of cores.")
    }
    
    # Optional: Check against the number of available cores
    available_cores <- parallel::detectCores()
    if (num_cores > available_cores) {
      warning(paste("Requested number of cores (", num_cores, ") exceeds the available cores (", available_cores, "). Proceed with caution."))
    }
    
    print(paste("Using", num_cores, "cores"))
  }
} else {
  num_cores <- 2
  print("Number of cores not provided, using default - 2 cores.")
}



###############################
#                             #
#   II.   F U N C T I O N S   #
#                             #
###############################

# Function to expand a single row
expand_row <- function(row, fasta_name_length) {
  genome_name_base <- row$genome_name
  copy_number <- row$copy_number
  genome_length <- fasta_name_length[[gsub("_.*", "", genome_name_base)]]
  
  # Create expanded data.table for a single row
  expanded_row_dt <- data.table(
    genome_name = paste0(genome_name_base, "_c", 1:copy_number),
    new_start = sample(1:genome_length, copy_number, replace = TRUE)
  )
  
  return(expanded_row_dt)
}

# Function to expand entire data table
expand_data_table_optimized <- function(dt, fasta_name_length, seed) {
  set.seed(seed)
  # Apply expand_row to each row and bind the results
  expanded_dt <- rbindlist(lapply(1:nrow(dt), function(i) expand_row(dt[i], fasta_name_length)))
  
  return(expanded_dt)
}

apply_binomial <- function(data, column_name, efficiency) {
  
  if (any(is.na(data[[column_name]]) | data[[column_name]] < 0 | data[[column_name]] != floor(data[[column_name]]))) {
    stop("Invalid value detected in column: ", column_name)
  }
  if (!is.numeric(efficiency) || efficiency < 0 || efficiency > 1) {
    stop("Efficiency must be a numeric value between 0 and 1.")
  }
  # Perform Bernoulli trial on the specified column
  data[[column_name]] <- rbinom(n = nrow(data), size = data[[column_name]], prob = efficiency)
  
  return(data)
}

#######################################################
#                                                     #
#   III.   F I N D I I N G  S A M P L E  T A B L E    #
#                                                     #
#######################################################


sample_table_path<-list.files(pattern = "_sample_table\\.csv$")
fasta_name_table<-read.csv("fasta_lengths.csv")
fasta_name_length <- split(fasta_name_table$fasta_lengths, fasta_name_table$fasta_names)


##########################################################################################
#                                                                                        #
#   IV.   S U B S E T T I N G  A  S A M P L E  F O R  F R A G/T A G M E N T A T I O N    #
#                                                                                        #
##########################################################################################

# Read each file as a data table and store in a list with appropriate names
data_tables <- setNames(lapply(sample_table_path, fread), gsub(pattern = ".csv$", replacement = "", x = sample_table_path))

# Apply the apply_binomial function to each data table
data_tables <- lapply(data_tables, apply_binomial, column_name = "copy_number", efficiency = as.numeric(fragmentation_fraction))


# Filter out rows where copy_number is 0
data_tables <- mclapply(data_tables, function(dt) {
  dt[dt$copy_number != 0, ]
}, mc.cores = num_cores)

# generate seeds for each file as they will be processed by different cores

if (exists("seed_value")) {
  set.seed(seed_value)
  worker_seeds <- sample(1:1e7,length(data_tables), replace = FALSE)
} else {
  #in the case there is no seed, each run will there be random seeds
  worker_seeds <- sample(1:1e7, length(data_tables), replace = FALSE)
}

#######################################################
#                                                     #
#   V.   C I R C U L A R ?  -  N E W  S T A R T S     #
#                                                     #
#######################################################

# making the new starting point of the genomes if it is circular, if it is not table is still made BUT the starting position is always 1!
if (circular_option == "Y") {
  print("Genome specifed as circular, each genome copy will get new starting position")
  # apply and get the table with new starta
  new_starts <- mclapply(seq_along(data_tables), function(i) expand_data_table_optimized(data_tables[[i]], fasta_name_length, worker_seeds[i]), mc.cores = num_cores)
  names(new_starts) <- names(data_tables)
} else {
  print("Genome specified as linear, all genome copies start with position 1")
  
  expand_rows_linear <- function(genome_name, copy_number) {
    expanded <- paste0(genome_name, "_c", 1:copy_number)
    data.table(genome_name = expanded, new_start = 1)
  }
  
  # Function to apply the expansion to a single table
  process_table_linear <- function(dt) {
    rbindlist(lapply(1:nrow(dt), function(i) {
      expand_rows_linear(dt$genome_name[i], dt$copy_number[i])
    }))
  }
  
  # Apply the process_table function to each table in the list
  new_starts <- lapply(data_tables, process_table_linear)
  
  
}

########################################################################
#                                                                      #
#   VI. E X T R A C T I N G  F R A G M E N T  C O O R D I N A T E S    #
#                                                                      #
########################################################################

# Sourcing the C++ script
sourceCpp('/usr/src/app/scripts/ampliseq/fragmentation.cpp')

# Use mclapply to run the loop in parallel
mclapply(seq_along(data_tables), function(i) {
  
  if (!is.null(worker_seeds)) {
    set.seed(worker_seeds[i])
  }
  
  # Get the table name
  table_name <- names(data_tables)[i]
  
  # Get the corresponding genome length
  prefix <- sub("^(.*?)_.*", "\\1", table_name)
  genome_length <- fasta_name_length[[paste0(prefix)]][1]
  
  # Call the C++ function
  expand_table_cpp(data_tables[[i]], table_name, genome_length, min(fragment_length), max(fragment_length), worker_seeds[i], output_folder)
}, mc.cores = num_cores)


# Function to save tables in a list to CSV files with modified file names
save_tables_to_csv <- function(tables_list, suffix, output_folder) {
  names(tables_list) <- gsub("_sample_table", suffix, names(tables_list)) # Rename the table names
  lapply(names(tables_list), function(name) {
    file_name <- file.path(output_folder, paste0(name, ".csv"))
    write.csv(tables_list[[name]], file_name, row.names = FALSE)
    cat("Saved table", name, "to file", file_name, "\n")
  })
}

# Save each table in data_tables with '_fragmentation_selected' suffix
save_tables_to_csv(data_tables, "_fragmentation_selected", output_folder)

# Save each table in new_starts with '_new_start' suffix
save_tables_to_csv(new_starts, "_new_start", output_folder)

print("Fragmentation is finished!")
print("In the fragmentation_selected.csv you will find the name of the genome copis that are fragmented")
print("In the new_start.csv you will find the start position of every genome copy included in fragmentation. If genome was linear, the start position for all genome copies is 1")



