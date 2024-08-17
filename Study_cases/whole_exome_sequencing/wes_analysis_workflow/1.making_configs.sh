#!/bin/bash

# Define arrays of parameter values
mutation_fraction=("0.4" "0.8")
num_reads=("35000000")
nr_copies=("1500" "3500") 

# Template configuration file path
template_config="config_iswes.yml"

# Loop through the parameter arrays and create a config file for each combination
for fraction in "${mutation_fraction[@]}"; do
    for reads in "${num_reads[@]}"; do
        for copies in "${nr_copies[@]}"; do
             # create a new config filename based on the parameters
             new_config="config_mut_fraction_${fraction}_r_${reads}_c_${copies}.yml"

             # copy the template config file to a new config file
             cp "$template_config" "$new_config"

             # replace parameters in the new config file
             sed -i "s/MUT_GENOME_FRACTION: \".*\"/MUT_GENOME_FRACTION: \"$fraction\"/" "$new_config"
             sed -i "s/NUM_READS: \".*\"/NUM_READS: \"$reads\"/" "$new_config"
             sed -i "s/NR_COPIES: \".*\"/NR_COPIES: \"$copies\"/" "$new_config"

             echo "Created $new_config"
        done
    done
done
