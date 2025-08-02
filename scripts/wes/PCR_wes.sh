#!/bin/bash


input_file=""
seed_value=-1 # Use -1 to indicate that the seed is not set
num_cores=2 # Default number of cores
polymerase_error_rate=1e-6 # Default polymerase error rate
num_cycles=25 # Default number of PCR cycles
midpoint_cycle=-1 # Calculated based on num_cycles if not directly provided
k_parameter=0.3 # Default k parameter for PCR
optimal_length_mode="default"
k_param_seq=0.05 # Default k parameter for sequencing
num_reads=-1 # Indicates that the number of reads is not provided

# Function to print usage
print_usage() {
  echo "Usage: $0 [options]"
  echo "--input_file {file} : Specify input file"
  echo "--set_seed {integer} : Sets a seed for reproducibility"
  echo "--cores {integer} : Defines number of cores"
  echo "--num_cycles {integer} : Define number of PCR cycles"
  echo "--midpoint_cycle {integer} : Define the midpoint cycle"
  echo "--k_parameter_pcr {float} : Define k parameter for PCR"
  echo "--optimal_length_mode {mode} : Define optimal length mode"
  echo "--k_parameter_seq {float} : Define k parameter for sequencing"
  echo "--num_reads {integer} : Define number of reads"
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --input_file)
      input_file="$2"
      shift # past argument
      shift # past value
      ;;
    --set_seed)
      seed_value="$2"
      shift # past argument
      shift # past value
      ;;
    --cores)
      num_cores="$2"
      shift # past argument
      shift # past value
      ;;
    --num_cycles)
      num_cycles="$2"
      shift # past argument
      shift # past value
      ;;
    --midpoint_cycle)
      midpoint_cycle="$2"
      shift # past argument
      shift # past value
      ;;
    --k_parameter_pcr)
      k_parameter="$2"
      shift # past argument
      shift # past value
      ;;
    --optimal_length_mode)
      optimal_length_mode="$2"
      shift # past argument
      shift # past value
      ;;
    --k_parameter_seq)
      k_param_seq="$2"
      shift # past argument
      shift # past value
      ;;
    --num_reads)
      num_reads="$2"
      shift # past argument
      shift # past value
      ;; 
    *)    # unknown option
      print_usage
      exit 1
      ;;
  esac
done


if [ -z "$input_file" ]; then
  echo "Input file not provided"
  exit 1
fi

if [ -n "$seed_value" ]; then
  echo "Seed set to $seed_value. "
else
  error "Seed not provided"
fi

echo "Using $num_cores cores"
echo "Polymerase error rate is set to $polymerase_error_rate"
echo "Defined number of PCR cycles: $num_cycles"
if [ "$midpoint_cycle" -ne -1 ]; then
  echo "Midpoint cycle defined as: $midpoint_cycle"
else
  midpoint_cycle=$((num_cycles * 60 / 100))
  echo "Midpoint cycle is not defined, using default where midpoint is a cycle corresponding to 60% of all cycles, which is now $midpoint_cycle"
fi
echo "k parameter defined as: $k_parameter"
echo "Optimal length mode: $optimal_length_mode"
if [ -n "$k_param_seq" ]; then
  echo "k parameter sequencing: $k_param_seq"
fi
if [ -n "$num_reads" ]; then
  echo "Number of reads: $num_reads"
else
  error "Number of reads not provided."
fi


# Extract the base path 
base_path=$(echo "$input_file" | sed -E 's|(.*)sample_.*|\1|')

# Extract sample_number and tech_rep_number using pattern matching
sample_number=$(echo "$input_file" | grep -o 'sample_[0-9]*' | head -1)
tech_rep_number=$(echo "$input_file" | grep -o 'tech_replicate_[0-9]*' | head -1)


# Extracting file name from the full file path
file_name=$(basename "$input_file")

echo "Currently processing $file_name"

# Construct the destination path
destination_path="${base_path}${sample_number}/${tech_rep_number}/"

# Construct the main PCR folder path
main_pcr_folder="${destination_path}PCR_reaction/"


# Call cpp script to run PCR
PCR_command="/usr/src/app/scripts/wes/PCR_cpp $input_file $num_cores $seed_value $num_cycles $midpoint_cycle $k_parameter $main_pcr_folder"

# Execute the command
eval $PCR_command

#Make the final table and abundance file names
final_table_name="${main_pcr_folder}PCR_all_amp_fragments.csv"
final_abundance_name="${main_pcr_folder}PCR_abundance.txt"

# Call the script to merge CSV files and generate the abundance file
/usr/src/app/scripts/wes/merging_csv_files "$main_pcr_folder" "$final_table_name" "$final_abundance_name"

# Check if optimal_length_mode is set to "NO"
if [ "$optimal_length_mode" != "NO" ]; then
    echo "making an additional file fasta_lengths.csv with the lengths of all amplicons"
    fasta_length_file_name="${main_pcr_folder}PCR_fasta_lengths.csv"
    
    /usr/src/app/scripts/wes/extracting_fasta_length.sh "$final_table_name" "$fasta_length_file_name"
fi

# Define the output file - read_pairs_per_fragment
output_nr_reads_per_fragments="${main_pcr_folder}read_pairs_per_fragment.csv"

# Base command
command_read_pairs="/usr/src/app/scripts/wes/calculate_read_pairs --abundance_file ${final_abundance_name} --optimal_length_mode \"${optimal_length_mode}\" --seed ${seed_value} --read_number ${num_reads} --output ${output_nr_reads_per_fragments}"

# Expand the command with k_value and fasta_lengths if optimal_length_mode is not "NO"
if [ "${optimal_length_mode}" != "NO" ]; then
    command_read_pairs+=" --k_value ${k_param_seq} --fasta_lengths ${fasta_length_file_name}"
fi

# Execute the command
eval "$command_read_pairs"

# Filter main pcr csv

filtered_amp_fragments_table_name="${main_pcr_folder}PCR_filtered_amp_fragments.csv"
command_filtering="/usr/src/app/scripts/wes/filtering_fragments_for_sequencing.sh $output_nr_reads_per_fragments $final_table_name $filtered_amp_fragments_table_name"

eval "$command_filtering"

# Sort the fragments by names

command_sorting="/usr/src/app/scripts/wes/sorting_frags_by_name.sh $filtered_amp_fragments_table_name"

eval "$command_sorting"

# Split table into number of cores

command_splitting="/usr/src/app/scripts/wes/splitting_filtered_fragments.sh $filtered_amp_fragments_table_name $num_cores"
eval "$command_splitting"

# preparing for processing fragments and making fasta files

temp_files=($(find "$main_pcr_folder" -maxdepth 1 -type f -name "worker_*_temp.csv" -print))
num_temp_files=${#temp_files[@]}

if [ ! -z "$seed_value" -a "$num_temp_files" -gt 0 ]; then
    # Call the R script, passing the seed_value and the number of seeds needed
    readarray -t worker_seeds_pcr < <(Rscript /usr/src/app/scripts/wes/generate_seed.R "$seed_value" "$num_temp_files")

    echo "Generated worker seeds:"
    for seed in "${worker_seeds_pcr[@]}"; do
        echo $seed
    done
else
    error "No seed value provided or no temp files to process; worker seeds will not be generated."
fi

process_file() {
    temp_file=$1
    worker_seed=$2
    sample_folder="${base_path}${sample_number}/"
    worker_id=$BASHPID
    temp_fasta_files="worker_${worker_id}_temp"
    temp_fasta_files_with_path="${main_pcr_folder}${temp_fasta_files}"
  
    /usr/src/app/scripts/wes/process_frags "$temp_file" "$sample_folder" "/usr/src/app/pipeline/SQL_database/" "$worker_seed" "$temp_fasta_files_with_path"
}
export main_pcr_folder base_path sample_number
export -f process_file

parallel --link -j "$num_cores" process_file ::: "${temp_files[@]}" ::: "${worker_seeds_pcr[@]}"

# Merge produced fasta files

# List temp FASTA files
temp_fasta_names=($(find "${main_pcr_folder}" -maxdepth 1 -type f -regextype posix-extended -regex ".*/worker_.*_temp-[0-9]+\.fasta"))

# Define the final FASTA file path
final_fasta="${main_pcr_folder}PCR_reaction.fasta"

echo "Merging temporary FASTA files into ${final_fasta}."

for temp_fasta_file in "${temp_fasta_names[@]}"; do
  cat "$temp_fasta_file" >> "$final_fasta"
  rm "$temp_fasta_file"
done

echo "Merging of FASTA files completed and temporary files removed."

# remerging the csv files and removal of the temp csv files

for temp_file in "${temp_files[@]}"; do
    cat "$temp_file" >> "${filtered_amp_fragments_table_name}"
done

for temp_file in "${temp_files[@]}"; do
    rm "$temp_file"
done

# clean up

rm "$final_abundance_name"

if [ "$optimal_length_mode" != "NO" ]; then
    rm $fasta_length_file_name
fi

# Prep for sequencing - split the file into the number of cores
/usr/src/app/scripts/wes/splitting_reads_csv.sh "$output_nr_reads_per_fragments" "$num_cores" "$main_pcr_folder"

# Prep for sequencing - split fasat file into the number of cores
/usr/src/app/scripts/wes/splitting_fasta_prep_seq.sh "$main_pcr_folder" "$final_fasta" "$main_pcr_folder"

csv_output_list="${main_pcr_folder}csv_list.txt"
fasta_output_list="${main_pcr_folder}fasta_list.txt"

# list the csv and fasta files - input for sequencing
find "${main_pcr_folder}" -maxdepth 1 -type f -name "temp_reads_*.csv" > "$csv_output_list"
find "${main_pcr_folder}" -maxdepth 1 -type f -name "temp_*.fasta" > "$fasta_output_list"

echo "CSV and FASTA file lists have been created."

# Merging temp-mutation-count.csv
temp_counts_names=($(find "${main_pcr_folder}" -maxdepth 1 -type f -regextype posix-extended -regex ".*/worker_.*_temp-mutation-count\.csv"))
final_count_name="${main_pcr_folder}seq_mutation_counts.csv"

# Merge all files into one
cat "${temp_counts_names[@]}" > "${final_count_name}.tmp"

# Process the merged file to sum duplicates
sort "${final_count_name}.tmp" | awk -F, '{
    if (last_key && last_key != $1 "," $2) {
        print last_key "," sum
        sum = 0
    }
    last_key = $1 "," $2
    sum += $3
}
END {
    if (last_key) print last_key "," sum
}' > "${final_count_name}"

# Remove temporary file
rm "${final_count_name}.tmp"

echo "Merging and summary complete. Results saved in ${final_count_name}"

for temp_counts_names in "${temp_counts_names[@]}"; do
    rm "$temp_counts_names"
done
