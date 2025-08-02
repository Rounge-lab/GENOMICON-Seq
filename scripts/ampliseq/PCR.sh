#!/bin/bash


# Initialize default values
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
frac_to_seq=-1 # Indicates that the fraction of PCR sample that goes to sequencing is not provided

# Function to print usage
print_usage() {
  echo "Usage: $0 [options]"
  echo "--input_file {file} : Specify input file"
  echo "--set_seed {integer} : Sets a seed for reproducibility"
  echo "--cores {integer} : Defines number of cores"
  echo "--polymerase_error_rate {float} : Define polymerase error rate"
  echo "--num_cycles {integer} : Define number of PCR cycles"
  echo "--midpoint_cycle {integer} : Define the midpoint cycle"
  echo "--k_parameter_pcr {float} : Define k parameter for PCR"
  echo "--optimal_length_mode {mode} : Define optimal length mode"
  echo "--k_parameter_seq {float} : Define k parameter for sequencing"
  echo "--num_reads {integer} : Define number of reads"
  echo "--frac_to_seq {float} : Define the fraction of PCR sample that will be sequenced"
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
    --polymerase_error_rate)
      polymerase_error_rate="$2"
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
    --frac_to_seq)
      frac_to_seq="$2"
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
echo "Optimal length mode: $optimal_length_mode."
if [ -n "$k_param_seq" ]; then
  echo "k parameter sequencing: $k_param_seq"
fi
if [ -n "$num_reads" ]; then
  echo "Number of reads: $num_reads. "
else
  error "Number of reads not provided."
fi


# Extract the base path up to "pipeline" and the next directory name, assuming "pipeline" is always part of the path
base_path=$(echo "$input_file" | sed -E 's|(.*)sample_.*|\1|')

# Extract sample_number and tech_rep_number using pattern matching
sample_number=$(echo "$input_file" | grep -o 'sample_[0-9]*' | head -1)
tech_rep_number=$(echo "$input_file" | grep -o 'tech_replicate_[0-9]*' | head -1)


# Extracting file name from the full file path
file_name=$(basename "$input_file")

# Execute the Python script with the full file path as an argument
python3 /usr/src/app/scripts/ampliseq/expanding_pcr_input_table.py "$input_file"

# Providing feedback to the user
echo "Currently processing $file_name."


# Construct the destination path for split files
destination_path="${base_path}${sample_number}/${tech_rep_number}/PCR/"

# Extract the PCR index number using sed with a regular expression
pcr_index_number=$(echo "$file_name" | sed -E 's/.*_([0-9]+)\.csv$/\1/')

# Construct the main PCR folder name using the extracted index
main_pcr_folder="${destination_path}PCR_${pcr_index_number}"


# Prepare the command to split the input file according to the number of cores
splitting_input_command="/usr/src/app/scripts/ampliseq/splitting_prep_pcr.sh $input_file $num_cores $destination_path"

# Execute the command
eval $splitting_input_command

# Extract the base file name by removing the '.csv' extension
file_name_base=$(echo "$file_name" | sed 's/\.csv$//')

# Correct approach to define the pattern without prepending the destination path
# Correct approach to define the pattern and ensure we correctly match and list the files
temp_files_base_pattern="${file_name_base}_[0-9]+\\.csv$"

# Use find to locate files in the destination path that match the pattern, ensuring the correct format for listing
list_temp_split_files=($(find "$destination_path" -maxdepth 1 -type f -regextype posix-extended -regex ".*${temp_files_base_pattern}" -printf "%p\n"))


# Ensure this line is correctly positioned in your script
num_temp_files=${#list_temp_split_files[@]}

if [ ! -z "$seed_value" -a "$num_temp_files" -gt 0 ]; then
    # Call the R script, passing the seed_value and the number of seeds needed
    readarray -t worker_seeds_pcr < <(Rscript /usr/src/app/scripts/ampliseq/generate_seed.R "$seed_value" "$num_temp_files")

    echo "Generated worker seeds:"
    for seed in "${worker_seeds_pcr[@]}"; do
        echo $seed
    done
else
    echo "No seed value provided or no temp files to process; worker seeds will not be generated."
fi


process_file() {
  temp_file=$1
  worker_seed=$2 
  temp_ampli_table_files="table_amplicons_${BASHPID}_temp"
  temp_ampli_table_with_path="${main_pcr_folder}/${temp_ampli_table_files}"

  /usr/src/app/scripts/ampliseq/PCR_cpp $temp_file $worker_seed $polymerase_error_rate $num_cycles $midpoint_cycle $k_parameter $frac_to_seq $temp_ampli_table_with_path
}

export -f process_file
export polymerase_error_rate num_cycles midpoint_cycle k_parameter main_pcr_folder frac_to_seq

# Use GNU parallel to run process_file on each file with its corresponding seed
parallel --link -j "$num_cores" process_file ::: "${list_temp_split_files[@]}" ::: "${worker_seeds_pcr[@]}"

rm "${list_temp_split_files[@]}"

# Construct the final table and abundance file names
final_table_name="${main_pcr_folder}/PCR_${pcr_index_number}_all_amp_fragments.csv"
final_abundance_name="${main_pcr_folder}/PCR_${pcr_index_number}_abundance.txt"

# Call the script to merge CSV files and generate the abundance file
/usr/src/app/scripts/ampliseq/merging_csv_files "$main_pcr_folder" "$final_table_name" "$final_abundance_name"


# Check if optimal_length_mode is set to "NO"
if [ "$optimal_length_mode" != "NO" ]; then
    echo "Making fasta_lengths.csv with the lengths of all amplicons"
    fasta_length_file_name="${main_pcr_folder}/PCR_${pcr_index_number}_fasta_lengths.csv"
    
    /usr/src/app/scripts/ampliseq/extracting_fasta_length.sh "$final_table_name" "$fasta_length_file_name"
fi


output_nr_reads_per_fragments="${main_pcr_folder}/PCR_${pcr_index_number}_reads_per_fragment.csv"

# Initialize the base command with proper quoting for the optimal_length_mode variable
command="/usr/src/app/scripts/ampliseq/calculate_read_pairs --abundance_file ${final_abundance_name} --optimal_length_mode \"${optimal_length_mode}\" --seed ${seed_value} --read_number ${num_reads} --output ${output_nr_reads_per_fragments}"

# Conditionally include k_value and fasta_lengths if optimal_length_mode is not "NO"
if [ "${optimal_length_mode}" != "NO" ]; then
    command+=" --k_value ${k_param_seq} --fasta_lengths ${fasta_length_file_name}"
fi

# Execute the command
eval "$command"

# Filter fragments that are going to be sequenced
filtered_amp_fragments_table_name="${main_pcr_folder}/PCR_${pcr_index_number}_filtered_amp_fragments.csv"

# Construct the command for filtering
command_filtering="/usr/src/app/scripts/ampliseq/preprocessing_pcr_csv.sh ${output_nr_reads_per_fragments} ${final_table_name} ${input_file} ${filtered_amp_fragments_table_name}"

# Execute the command for filtering
eval "$command_filtering"

echo "Filtering completed. Filtered fragments table is at: ${filtered_amp_fragments_table_name}"

# sorting the filtered fragments

sort -o "$filtered_amp_fragments_table_name" "$filtered_amp_fragments_table_name"

# Split table into number of cores

command_splitting="/usr/src/app/scripts/ampliseq/splitting_filtered_fragments.sh $filtered_amp_fragments_table_name $num_cores"
eval "$command_splitting"

# preparing for processing fragments and making fasta files

temp_files=($(find "$main_pcr_folder" -maxdepth 1 -type f -name "worker_*_temp.csv" -print))
num_temp_files=${#temp_files[@]}

if [ ! -z "$seed_value" -a "$num_temp_files" -gt 0 ]; then
    # Call the R script, passing the seed_value and the number of seeds needed
    readarray -t worker_seeds_fasta < <(Rscript /usr/src/app/scripts/ampliseq/generate_seed.R "$seed_value" "$num_temp_files")

    echo "Generated worker seeds:"
    for seed in "${worker_seeds_fasta[@]}"; do
        echo $seed
    done
else
    error "No seed value provided or no temp files to process; worker seeds will not be generated."
fi

process_temp_file() {
  temp_file_csv=$1
  worker_seed=$2
  sample_folder="${base_path}${sample_number}/"
  SQL_database_folder="/usr/src/app/pipeline/SQL_database/"
  length_file="${base_path}${sample_number}/fasta_lengths.csv"
  worker_id=$BASHPID
  temp_fasta_file="worker_${worker_id}_temp"
  temp_fasta_files_with_path="${main_pcr_folder}/${temp_fasta_file}"

  /usr/src/app/scripts/ampliseq/process_fragment_is_ampliseq "$temp_file_csv" "$sample_folder" "$SQL_database_folder" "$worker_seed" "$length_file" "$temp_fasta_files_with_path"

} 
export base_path sample_number main_pcr_folder
export -f process_temp_file

parallel --link -j "$num_cores" process_temp_file ::: "${temp_files[@]}" ::: "${worker_seeds_fasta[@]}"

# List temp FASTA files
temp_fasta_names=($(find "${main_pcr_folder}" -maxdepth 1 -type f -regextype posix-extended -regex ".*/worker_.*_temp-[0-9]+\.fasta"))

# Define the final FASTA file path
final_fasta="${main_pcr_folder}/PCR_${pcr_index_number}.fasta"

echo "Merging temporary FASTA files into ${final_fasta}."

for temp_fasta_file in "${temp_fasta_names[@]}"; do
  cat "$temp_fasta_file" >> "$final_fasta"
  rm "$temp_fasta_file"
done

echo "Merging of FASTA files completed and temporary fasta files removed."

# clean up fasta_length_file_name final_abundance_name and final_table_name

rm "$final_abundance_name"

rm "$final_table_name"

if [ "$optimal_length_mode" != "NO" ]; then
    rm $fasta_length_file_name
fi

# merge the temp_files (split filtered_amp_fragments_table_name) and remake it
for temp_file in "${temp_files[@]}"; do
    cat "$temp_file" >> "${filtered_amp_fragments_table_name}"
done

for temp_file in "${temp_files[@]}"; do
    rm "$temp_file"
done

echo "Merging PCR mutations tracking temporary csv files into a final csv file"

# Define the output csv file name
pcr_mutations="${main_pcr_folder}/PCR_${pcr_index_number}-polymerase_error_mutations.csv"

# List temp csv files matching the pattern
temp_csv_names=($(find "${main_pcr_folder}" -maxdepth 1 -type f -regextype posix-extended -regex ".*/worker_[0-9]+_temp-mutation_tracking\.csv"))

# Initialize a flag to check if the header has been copied
header_copied=false

# Check if temp csv names array is not empty
if [ ${#temp_csv_names[@]} -gt 0 ]; then
    # Loop through each file found
    for file in "${temp_csv_names[@]}"; do
        # Check if the header has been copied
        if [ "$header_copied" = false ]; then
            # Copy the first file with header
            cat "$file" > "$pcr_mutations"
            # Set the flag indicating the header has been copied
            header_copied=true
        else
            # Append the file without the header
            tail -n +2 "$file" >> "$pcr_mutations"
        fi
    done
    echo "Merge complete. Counting the occurance of the unique mutations."

    # Consolidate records by summing copy_number for identical pcr_mutations
    awk -F, 'NR == 1 { print $0; next } { count[$1] += $2 } END { for (mutation in count) print mutation, count[mutation] }' OFS=, "$pcr_mutations" > "${pcr_mutations}.tmp" && mv "${pcr_mutations}.tmp" "$pcr_mutations"

    echo "Counting complete. Removing temporary files..."
    
    # Remove the original temporary files
    rm "${temp_csv_names[@]}"
    
    echo "Temporary files removed."
else
    # No files found, so echo the message and create an empty log file
    echo "No PCR mutations have been introduced."
    touch "${main_pcr_folder}/PCR_${pcr_index_number}_no_pcr_mutations.log"
fi


# Merge temp files with inserted real mutations into a final csv file

echo "Merging true inserted mutation tracking files into a final file"

# Define the output csv file name
real_mutations="${main_pcr_folder}/PCR_${pcr_index_number}-real_inserted_mutations.csv"

# List temp csv files matching the pattern
temp_csv_names_real_mut=($(find "${main_pcr_folder}" -maxdepth 1 -type f -regextype posix-extended -regex ".*/worker_[0-9]+_temp-real_mutations\.csv"))

# Initialize a flag to check if the header has been copied
header_copied=false

# Check if temp csv names array is not empty
if [ ${#temp_csv_names_real_mut[@]} -gt 0 ]; then
    # Loop through each file found
    for file in "${temp_csv_names_real_mut[@]}"; do
        # Check if the header has been copied
        if [ "$header_copied" = false ]; then
            # Copy the first file with header
            cat "$file" > "$real_mutations"
            # Set the flag indicating the header has been copied
            header_copied=true
        else
            # Append the file without the header
            tail -n +2 "$file" >> "$real_mutations"
        fi
    done
    echo "Merge complete. Coutning the occurance of the unique mutations."

    # Consolidate records by summing copy_number for identical pcr_mutations
    awk -F, 'NR == 1 { print $0; next } { count[$1 "," $2] += $3 } END { for (key in count) print key, count[key] }' OFS=, "$real_mutations" > "${real_mutations}.tmp" && mv "${real_mutations}.tmp" "$real_mutations"

    echo "Counting complete. Removing temporary files..."
    
    # Remove the original temporary files
    rm "${temp_csv_names_real_mut[@]}"
    
    echo "Temporary files removed."
else
    # No files found, so echo the message and create an empty log file
    echo "No true mutations have been introduced."
    touch "${main_pcr_folder}/PCR_${pcr_index_number}_no_introduced_true_mutations.log"
fi



echo "Preparing input for the sequencing."

# Prep for sequencing - split the file into the number of cores
/usr/src/app/scripts/ampliseq/splitting_reads_csv.sh "$output_nr_reads_per_fragments" "$num_cores" "$main_pcr_folder"

# Split FASTA
/usr/src/app/scripts/ampliseq/splitting_fasta_prep_seq.sh "$main_pcr_folder" "$final_fasta" "$main_pcr_folder"

csv_read_files="${main_pcr_folder}/temp_reads_*.csv"
csv_output_list="${main_pcr_folder}/PCR_${pcr_index_number}_csv_list.txt"
fasta_files="${main_pcr_folder}/temp_*.fasta"
fasta_output_list="${main_pcr_folder}/PCR_${pcr_index_number}_fasta_list.txt"

# Use find command to correctly handle paths and generate lists
find "${main_pcr_folder}" -maxdepth 1 -type f -name "temp_reads_*.csv" > "$csv_output_list"
find "${main_pcr_folder}" -maxdepth 1 -type f -name "temp_*.fasta" > "$fasta_output_list"

echo "CSV and FASTA file lists have been created."