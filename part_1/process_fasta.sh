#!/bin/bash

# define the input fasta file
input_fasta="/Path/to/modified_mature_sequences.fas"

# directory to store output files
output_dir="mirna_family_single_fasta"

# create directory (if not existent)
mkdir -p "$output_dir"

# loop over each CSV file in the current directory
for csv_file in *.csv; do
    # get the base name without extension
    base_name=$(basename "$csv_file" .csv)
    # assign the output file name
    output_file="$output_dir/${base_name}.fas"
    # use seqkit grep to filter sequences and save to the output file
    seqkit grep -f "$csv_file" "$input_fasta" > "$output_file"
done

echo "All done! Files are saved in $output_dir"
