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

###############################################################################################################
#                                                                                                             #
#      0.   W O R K I N G  D I R E C T O R Y  & S E T T I N G  A  S E E D  &  N U M B E R  O F  C O R E S     #
#                                                                                                             #
###############################################################################################################

if (any(grepl("--output", commandArgs()))) { 
  output_arg <- which(grepl("--output", commandArgs())) 
  if (length(output_arg) > 0 && output_arg < length(commandArgs())) { 
    output_folder <- commandArgs()[output_arg + 1] 
    print(paste("Ok, we will use this", output_folder, "as the output folder, ok!")) } } 

setwd(output_folder) 

# argument --set_seed {integer} - sets a seed for reproducibility of random processes 
if (any(grepl("--set_seed", commandArgs()))) {
  seed_arg <- which(grepl("--set_seed", commandArgs()))
  if (length(seed_arg) > 0 && seed_arg < length(commandArgs())) {
    seed_value <- as.integer(commandArgs()[seed_arg + 1])
    if (!is.na(seed_value) && seed_value >= 0) {
      set.seed(seed_value)
      print(paste0("Specified seed value: ", seed_value,""))
    } else {
      stop("Invalid seed value, must be an integer")
    }
  }
} else {
  print("No seed were provided")
}

# argument --cores {integer} defines number of cores, if not given process is run on 1 core sequentially
if (any(grepl("--cores", commandArgs()))) {
  core_arg <- which(grepl("--cores", commandArgs()))
  if (length(core_arg) > 0 && core_arg < length(commandArgs())) {
    num_cores <- as.integer(commandArgs()[core_arg + 1])
    
    if (is.na(num_cores) || num_cores <= 0 || num_cores != floor(num_cores)) {
      stop("Invalid number of cores.")
    }
    
    available_cores <- parallel::detectCores()
    if (num_cores > available_cores) {
      warning(paste("Requested number of cores (", num_cores, ") exceeds the available cores (", available_cores, ")"))
    }
    
    print(paste("Specified number of cores:", num_cores))
  }
} else {
  num_cores <- 2
  print("Number of cores not provided, using default - 2 cores")
}


##################################################
#                                                #
#      I.   D E F I N I N G   S E Q U E N C E    #
#                                                #
##################################################

# Load the fasta sequence
if (any(grepl("^--fasta_sequence$", commandArgs()))) {
  # If fasta_sequence argument is provided, extract the file name
  fasta_arg_index <- which(grepl("^--fasta_sequence$", commandArgs(trailingOnly = TRUE)))
  
  if (fasta_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    fasta_file <- commandArgs(trailingOnly = TRUE)[fasta_arg_index + 1]
    fasta_path <- file.path("/usr/src/app/pipeline/", fasta_file)
    # Check if the file exists
    if (file.exists(fasta_path)) {
      # Read the file
      fasta_content <- readDNAStringSet(fasta_path)
      print("Fasta file loaded")
    } else {
      stop("The specified fasta file does not exist")
    }
  }
} else {
  stop("Fasta sequence needs to be specified")
}

fasta_names <- names(fasta_content)
fasta_lengths <- nchar(fasta_content)
fasta_name_length_table<-data.table(fasta_names, fasta_lengths)
write.csv(fasta_name_length_table, "fasta_lengths.csv", row.names = F)
print("fasta_length.csv - generated.")

################################################################################################
#                                                                                              #
#      II. N U C L E O T I D E  C O N T E X T  &  M U T A T I O N  S I G N A T U R E           #
#                                                                                              #
################################################################################################
if (any(grepl("--mutational_context", commandArgs()))) {
  
  mut_context_pos <- which(grepl("--mutational_context", commandArgs()))
  if (length(mut_context_pos) > 0 && mut_context_pos < length(commandArgs())) {
    mut_context <- unlist(strsplit(commandArgs()[mut_context_pos + 1], ","))
    
    char_lengths <- nchar(mut_context)
    all_same_length <- all(char_lengths == char_lengths[1])
    valid_lengths <- all(char_lengths %in% c(1, 3, 5))
    
    if (all_same_length && valid_lengths) {
      print(paste("Specified target for mutational process: ",mut_context))
    } else {
      stop("Invalid mutational context. All elements must have the same length, and the length must be 1, 3, or 5 characters.")
    }
  }
}


#Define the Transitions vs Transversions ratio --Ts_Tv_ratio {Ts:Tv} 
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
        
        print(paste0("Specified transition vs Transversion ratio: ", ts_tv_ratio))
      } else {
        stop("Invalid format for --Ts_Tv_ratio. Expected format 'integer:integer'.")
      }
    } else {
      stop("Invalid format for --Ts_Tv_ratio. Expected format 'integer:integer'.")
    }
  }
} else if (any(grepl("^--variant_probability_table$", commandArgs()))) {
  #load variant probability table if specified
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
  print("table with variant mutation probability not provided, each nucleotide can mutate to other 3 variants with equal probability")
  variant_prob_table<-data.frame(
    Nucleotide=c("A","A","A","T","T","T","C","C","C","G","G","G"),
    Variant=c("T","C","G","A","C","G","G","A","T","C","A","T"),
    Probability=rep(100/3)
  )
}




##############################################################
#                                                            #
#         III. W H A T  I S  I N  T H E  S A M P L E ?       #
#                                                            #
##############################################################


#If multifasta file is provided by default the copy number applies to all the fasta files 
if (any(grepl("--nr_copies", commandArgs()))) {
  arg_pos <- which(grepl("--nr_copies", commandArgs()))
  
  if (length(arg_pos) > 0 && arg_pos < length(commandArgs())) {
    nr_copies <- as.integer(commandArgs()[arg_pos + 1])
    
    if (is.na(nr_copies) || nr_copies <= 0 || nr_copies != as.integer(nr_copies)) {
      stop("Invalid number of copies. It must be a positive integer")
    }
    
    print(paste0("Number of copies set to ", nr_copies))
    nr_copies <- rep(nr_copies, times = length(fasta_content))
  } 
} else if(any(grepl("^--multifasta_copies$", commandArgs()))) {
  multifasta_copies_table_arg_index <- which(grepl("^--multifasta_copies$", commandArgs(trailingOnly = TRUE)))
  if (multifasta_copies_table_arg_index < length(commandArgs(trailingOnly = TRUE))) {
    multifasta_copies_table_file <- commandArgs(trailingOnly = TRUE)[multifasta_copies_table_arg_index + 1]
    multifasta_copies_table_path <- file.path("/usr/src/app/pipeline/", multifasta_copies_table_file)
    if (file.exists(multifasta_copies_table_path)) {
      multifasta_copies_table<-read.csv(multifasta_copies_table_path)
      nr_copies<-as.numeric(multifasta_copies_table$copy_number)
      print("Table with a copy number for each multifasta file is loaded")
    }else {
      stop("The specified table not found.")
    } 
  } 
  
} else {
  print("No copy number provided, default - 1000 copies, if multifasta file provided, 1000 copies of each fasta sequence will be used.")
  nr_copies <- 1000
  nr_copies<-rep(nr_copies, time=length(fasta_content))
}


#################################################################################################
#                                                                                               #
#         IV.  D E F I N I N G  N E C E S S A R Y  F U N C T I O N S  &  V A R I A B L E S      #
#                                                                                               #
#################################################################################################
#Count the length of the mutational context
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
    stop("Length of mutational context did not met the length requirement, 1, 3 or 5, operating with either single nuclotides, tri or pentacontexts!")
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

#Function for finding the rows that matche the context
querying_SQL_with_context <- function(fasta_header, SQL_column, SQL_column_rev, mut_context) {
  sql_context <- paste0("'", mut_context, "'", collapse = ", ")
  
  conn <- dbConnect(SQLite(), dbname = paste0("/usr/src/app/pipeline/SQL_database/",fasta_header, ".sqlite"))
  
  query_5 <- sprintf(paste0("SELECT ROWID as row_number FROM ", fasta_header, " WHERE ", SQL_column, " IN (%s)"), sql_context)
  result_5 <- dbGetQuery(conn, query_5)
  mut_positions_5 <- result_5$row_number
  
  query_3 <- sprintf(paste0("SELECT ROWID as row_number FROM ", fasta_header, " WHERE ", SQL_column_rev, " IN (%s)"), sql_context)
  result_3 <- dbGetQuery(conn, query_3)
  mut_positions_3 <- result_3$row_number
  
  dbDisconnect(conn)
  
  return(list(mut_positions_5 = mut_positions_5, mut_positions_3 = mut_positions_3))
}
#Function that subselects the SQL database table, adds column variant and mutates it to desired variant
get_mutations_table <- function(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants) {
  conn <- dbConnect(SQLite(), dbname = paste0("/usr/src/app/pipeline/SQL_database/",fasta_header, ".sqlite"))
  
  row_ids_str_5 <- paste(mut_positions_5_selected, collapse = ", ")
  query_5 <- sprintf(paste0("SELECT * FROM ", fasta_header, " WHERE ROWID IN (%s)"), row_ids_str_5)
  mutations_table_5 <- dbGetQuery(conn, query_5)
  mutations_table_5$position<-mut_positions_5_selected
  mutations_table_5$Variant <- sapply(mutations_table_5$nucleotide, generate_variants)
  
  row_ids_str_3 <- paste(mut_positions_3_selected, collapse = ", ")
  query_3 <- sprintf(paste0("SELECT * FROM ", fasta_header, " WHERE ROWID IN (%s)"), row_ids_str_3)
  mutations_table_3 <- dbGetQuery(conn, query_3)
  mutations_table_3$position<-mut_positions_3_selected
  mutations_table_3$Variant <- sapply(mutations_table_3$comp_nucleotide, generate_variants)
  mutations_table_3$Variant <- chartr("ACGT", "TGCA", mutations_table_3$Variant)
  
  dbDisconnect(conn)
  
  return(rbind(mutations_table_5, mutations_table_3))
}                       

### Function that inserts the mutations from a sample_table into the sequences
mutate_sequence <- function(seq, mutation_string) {
  if (mutation_string == "No mutation") {
    return(seq)
  }
  
  mutation_list <- strsplit(mutation_string, ",")[[1]]
  positions <- as.numeric(sapply(mutation_list, function(x) unlist(strsplit(x, "[A-Z]"))[1]))
  nucleotides <- sapply(mutation_list, function(x) substr(x, nchar(x), nchar(x)))
  
  if(any(positions > length(seq)) | any(positions < 1)){
    stop(paste("Invalid positions detected for mutation string:", mutation_string))
  }
  
  mutated_seq <- replaceAt(seq, at = IRanges(start=positions, end=positions), value=DNAStringSet(nucleotides))
  
  return(mutated_seq)
}

#######################################################################################################
#                                                                                                     #
#         V.  M U T A T E D  G E N O M E   F R A C T I O N   O R   M U T A T I O N   R A T E ?        #
#                                                                                                     #
#######################################################################################################

###this will trigger either mutation rate or mutation frequency part of the script, one of those has to be given
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
    print("Fraction of mutation not specified, using default of 0.05 (5%)")
    pos_fraction <- 0.05
  }
  
} else if (any(grepl("--specific_mut_rate", commandArgs()))) {
  print("Fraction of mutated genome not provided, using specified mutation rate")
  arg_pos <- which(grepl("--specific_mut_rate", commandArgs()))
  if (length(arg_pos) > 0 && arg_pos < length(commandArgs())) {
    specific_mut_rate <- as.numeric(commandArgs()[arg_pos + 1])
    if (is.na(specific_mut_rate) || specific_mut_rate < 0 || specific_mut_rate > 1) {
      stop("Invalid specific mutation rate")
    }
    print(paste0("Specific mutation rate:", specific_mut_rate))
  } else {
    print("Specific mutation rate not defined, continuing with mutated genome fraction.")
    specific_mut_rate <- NULL
  }
  
} else {
  stop("None of the genome mutation options was provided, at least one of them has to be specified!")
}

################################################################################
#                                                                              #
#           VA    M U T A T E D   G E N O M E   F R A C T I O N                #  
#                                                                              #
################################################################################

##run in parallel for each sequence in the fasta_content
cl<-makeCluster(num_cores)
registerDoParallel(cl)  

worker_seeds <- NULL

if (exists("seed_value")) {
  # Generate worker seeds since main seed is set
  count_content <- length(fasta_content)
  worker_seeds <- sample(1:1e9, count_content, replace = FALSE)
} 


foreach(seq_index = 1:length(fasta_content), .packages = c('dplyr', 'Biostrings', 'DBI', 'RSQLite', 'data.table')) %dopar% {
  
  fasta_header <- names(fasta_content)[seq_index]

    if (!is.null(worker_seeds)) {
    worker_seed <- worker_seeds[seq_index]
    set.seed(worker_seed)
  }
  
    if (exists("mut_genome")) {
    if (exists("mut_context")) {

      mut_positions_3_5 <- querying_SQL_with_context(fasta_header, SQL_column, SQL_column_rev, mut_context)
      mut_positions_5 <- mut_positions_3_5$mut_positions_5
      mut_positions_3 <- mut_positions_3_5$mut_positions_3
      mut_positions <- sort(unique(c(mut_positions_5, mut_positions_3)))      
    } else {
      mut_positions<-1:width(fasta_content[seq_index])
    }
    
    tot_mut_number<-round(length(mut_positions)*pos_fraction)
    

    picked_mut_positions <- sample(mut_positions, tot_mut_number,replace = FALSE)
    picked_mut_positions <- sort(picked_mut_positions)
    
    
    if (exists("mut_positions_5")) {
      
      
      mut_positions_5_selected<-picked_mut_positions[picked_mut_positions %in% mut_positions_5]
      mut_positions_3_selected<-picked_mut_positions[picked_mut_positions %in% mut_positions_3]
      
      mutations_table <- get_mutations_table(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants)
      
    } else { 
      mut_positions_5_selected<-sample(picked_mut_positions, size = round(0.5 * length(picked_mut_positions)), replace = FALSE)
      mut_positions_3_selected <- setdiff(picked_mut_positions, mut_positions_5_selected)
      
      mutations_table <- get_mutations_table(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants)
    }
       
    #NMGs - non-mutated genomes
    MGs_number<-round(nr_copies[seq_index]*mut_genome)
    
    sample_table_NMG<-data.frame(
      genome_name=paste0(fasta_header, "_NMG"),
      variant_position = "No mutation",
      copy_number = nr_copies[seq_index] - MGs_number 
    )
    
    
    #MGs - mutated genomes
    sample_table_MG<-data.frame(
      genome_name=rep("",times = MGs_number),
      variant_position=rep("",times = MGs_number),
      copy_number = rep("",times = MGs_number)
    )
    
    mutations_table$pos_var<-paste0(mutations_table$position, mutations_table$Variant)
    all_possible_mutations<-mutations_table$pos_var
    
    for (i in 1:MGs_number) {
      num_mutations <- sample(length(all_possible_mutations), 1)
      
      selected_mutations <- sample(all_possible_mutations, num_mutations, replace = FALSE)
      selected_mutations <- selected_mutations[order(as.numeric(gsub("[A-Za-z]", "", selected_mutations)))]
      selected_mutations <- paste(selected_mutations, collapse = ",")
      
      
      sample_table_MG$genome_name[i]<-fasta_header
      sample_table_MG$variant_position[i]<-selected_mutations
      sample_table_MG$copy_number[i]<-1
      
    }
    
    sample_table_MG <- sample_table_MG %>%
      group_by(variant_position) %>%
      summarise(genome_name = dplyr::first(genome_name), copy_number = n()) %>%
      ungroup() %>% mutate(genome_name=paste0(genome_name, "_MG",row_number())) %>% 
      select(genome_name, variant_position, copy_number)
    
    
    sample_table<-rbind(sample_table_NMG, sample_table_MG)
    

    write.csv(mutations_table, paste0(fasta_header, "_mutations_table.csv"), row.names = F)
    write.csv(sample_table, paste0(fasta_header, "_sample_table.csv"), row.names = F)
    
    
  }else {
    fasta_header <- names(fasta_content)[[seq_index]]
    
    ################################################################################
    #                                                                              #
    #           VB   S P E C I F I C   M U T A T I O N   R A T E                   #     
    #                                                                              #
    ################################################################################
    
    if (exists("specific_mut_rate")) {

      fasta_length<-length(fasta_content[[seq_index]])
      lambda<-specific_mut_rate * fasta_length
      
      prob_MGs<-rpois(n=nr_copies[seq_index], lambda = lambda)
      
      has_mutation<-prob_MGs > 0
      n_mutations<-prob_MGs[has_mutation]
      
      if(sum(n_mutations) > 0){
        
        MGs_number<-length(n_mutations)
        sample_table_NMG<-data.frame(
          genome_name=paste0(fasta_header, "_NMG"),
          variant_position = "No mutation",
          copy_number = nr_copies[seq_index] - MGs_number 
        )

        MGs_table<-data.frame(
          genome_name = paste0(fasta_header, "_MG", 1:MGs_number),
          nr_mutations = n_mutations,
          positions = rep("", times = MGs_number),
          variant_position = rep("", times = MGs_number)
        )
        if (exists("mut_context")) {
          mut_positions_3_5 <- querying_SQL_with_context(fasta_header, SQL_column, SQL_column_rev, mut_context)
          mut_positions_5 <- mut_positions_3_5$mut_positions_5
          mut_positions_3 <- mut_positions_3_5$mut_positions_3
          mut_positions <- sort(unique(c(mut_positions_5, mut_positions_3)))

          MGs_table$positions <- sapply(MGs_table$nr_mutations, function(x) {
            paste(sample(mut_positions, x, replace = F), collapse = ",")
          })
          picked_mut_positions<-unique(sort(as.numeric(unlist(strsplit(MGs_table$positions, ",")))), decreasing = FALSE)

          mut_positions_5_selected<-picked_mut_positions[picked_mut_positions %in% mut_positions_5]
          mut_positions_3_selected<-picked_mut_positions[picked_mut_positions %in% mut_positions_3]
          
          mutations_table <- get_mutations_table(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants)
          
          
        } else {
          mut_positions<-1:length(fasta_content[[seq_index]])
          MGs_table$positions <- sapply(MGs_table$nr_mutations, function(x) {
            paste(sample(mut_positions, x, replace = F), collapse = ",")
          })
          picked_mut_positions<-paste(unique(sort(as.numeric(unlist(strsplit(MGs_table$positions, ",")))), decreasing = FALSE), collapse = ",")
          
          mut_positions_5_selected<-sample(picked_mut_positions, size = round(0.5 * length(picked_mut_positions)), replace = FALSE)
          mut_positions_3_selected <- as.numeric(unlist(strsplit(setdiff(picked_mut_positions, mut_positions_5_selected), ",")))
          
          mutations_table <- get_mutations_table(fasta_header, mut_positions_5_selected, mut_positions_3_selected, generate_variants)
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
        
        write.csv(mutations_table, paste0(fasta_header, "_mutations_table.csv"), row.names = F)
        write.csv(sample_table, paste0(fasta_header, "_sample_table.csv"), row.names = F)
        
      } else {
        print(paste("No copies of genome",fasta_header,"were mutated. The sample table with only non-mutated genomes will be made for this genome, copy number will correspond to the number of copies"))
        sample_table<-data.frame(
          genome_name=paste0(fasta_header, "_NMG"),
          variant_position = "No mutation",
          copy_number = nr_copies[seq_index]
        )
        write.csv(sample_table, paste0(fasta_header, "_sample_table.csv"), row.names = F)
        
      }
    }
  }
}

print("Writting out the outputs!")


stopImplicitCluster()

print("Process finished")
print("*_mutations_table.csv contains the information about mutations that was insreted into mutated genoms (one file per sequence in fasta file)")
print("*_sample_table.csv contains the information about how many genomes got a mutation and which positions/nucleotide were mutated")
print("Be aware, if the number of copies is too low as well as the mutation rate (if it was chosen, if not, then ignore this warning), it might happen that there are no mutated genomes at all, check and adjust the parameters!")
