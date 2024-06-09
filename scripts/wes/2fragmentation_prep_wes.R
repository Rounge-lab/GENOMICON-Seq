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
    print(paste("Using", input_folder, "as the input folder")) 
  } 
} 
setwd(input_folder)

if (any(grepl("--output_directory", commandArgs()))) { 
  output_arg <- which(grepl("--output_directory", commandArgs())) 
  if (length(output_arg) > 0 && output_arg < length(commandArgs())) { 
    output_folder <- commandArgs()[output_arg + 1] 
    print(paste("Using", output_folder, "as the output folder")) 
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


##Define proportion of a sample in % that will go to fragmentation/tagmentation
if (any(grepl("--fragmentation_fraction", commandArgs()))) {
  fragmentation_fraction_pos <- which(grepl("--fragmentation_fraction", commandArgs()))
  
  if (length(fragmentation_fraction_pos) > 0 && fragmentation_fraction_pos < length(commandArgs())) {
    fragmentation_fraction <- as.numeric(commandArgs()[fragmentation_fraction_pos + 1])
    
    # Check if fragmentation_fraction is between 0 and 1
    if (is.na(fragmentation_fraction) || fragmentation_fraction < 0 || fragmentation_fraction > 1) {
      stop("Invalid fragmentation fraction. The value must be between 0 and 1.")
    }
    
    print(paste0("Fragmentation fraction: ", fragmentation_fraction))
  }
} else {
  print("Fragmentation fraction not defined, using default - 5%")
  fragmentation_fraction <-1
}






###############################
#                             #
#   II.   F U N C T I O N S   #
#                             #
###############################


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



###################################################################################
#                                                                                 #
#   IV.   S U B S E T T I N G  A  S A M P L E  F O R  F R A GM E N T A T I O N    #
#                                                                                 #
###################################################################################

# Read each file as a data table and store in a list with appropriate names
data_tables <- setNames(lapply(sample_table_path, fread), gsub(pattern = ".csv$", replacement = "", x = sample_table_path))

# Apply the apply_binomial function to each data table
data_tables <- lapply(data_tables, apply_binomial, column_name = "copy_number", efficiency = as.numeric(fragmentation_fraction))

# Filter out rows where copy_number is 0
data_tables <- lapply(data_tables, function(dt) {
  dt[dt$copy_number != 0, ]
})



########################################################################
#                                                                      #
#   VI. E X T R A C T I N G  F R A G M E N T  C O O R D I N A T E S    #
#                                                                      #
########################################################################


# Function to save tables in a list to CSV files with modified file names
save_tables_to_csv <- function(tables_list, suffix, output_folder) {
  names(tables_list) <- gsub("_sample_table", suffix, names(tables_list)) # Rename the table names
  lapply(names(tables_list), function(name) {
    file_name <- file.path(output_folder, paste0(name, ".csv"))
    fwrite(tables_list[[name]], file_name, sep = ";", row.names = FALSE)
    cat("Saved table", name, "to file", file_name, "\n")
  })
}

# Save each table in data_tables with '_fragmentation_selected' suffix
save_tables_to_csv(data_tables, "_fragmentation_selected", output_folder)




