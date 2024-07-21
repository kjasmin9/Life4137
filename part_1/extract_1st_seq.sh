#!/bin/bash

# create empty output file
output_file="all_fasta_1seq.fas" > $output_file

# loop through each fasta file in the directory adn extract first ID and sequence
for file in *.fas
do
  # output in output_file
  head -n 2 "$file" >> $output_file
done
