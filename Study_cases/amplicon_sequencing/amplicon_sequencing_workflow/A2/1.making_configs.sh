#!/bin/bash

# Define arrays of parameter values
mutation_rates=("0" "5e-6" "5e-5" "5e-4")
num_reads=("1500000")
nr_copies=("200000")
poly_error_rate=("0" "1e-6" "1e-5" "1e-4")

# Template config file path
template_config="config_ampliseq.yml"

for rate in "${mutation_rates[@]}"; do
    for reads in "${num_reads[@]}"; do
        for copies in "${nr_copies[@]}"; do
	        for poly_rate in "${poly_error_rate[@]}"; do
                # create new config filename based on the parameters
                new_config="config_mut_r_${rate}_r_${reads}_c_${copies}_poly_r_${poly_rate}.yml"

                # copy template config file to a new config file
                cp "$template_config" "$new_config"

                # replace parameters in the new config file
                sed -i "s/MUTATION_RATE: \".*\"/MUTATION_RATE: \"$rate\"/" "$new_config"
                sed -i "s/NUM_READS: \".*\"/NUM_READS: \"$reads\"/" "$new_config"
                sed -i "s/NR_COPIES: \".*\"/NR_COPIES: \"$copies\"/" "$new_config"
	            sed -i "s/POL_ERROR_RATE: \".*\"/POL_ERROR_RATE: \"$poly_rate\"/" "$new_config"

                echo "Created $new_config"
             done
        done
    done
done