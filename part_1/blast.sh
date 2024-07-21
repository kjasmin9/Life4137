#!/bin/bash
#SBATCH --partition=defq 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50g
#SBATCH --time=24:00:00
#SBATCH --job-name=blastn_for_430_fas_files

#load module
module load blast-uoneasy/2.14.1-gompi-2023a

# define input and output directories
input_dir="4blast_fas"
output_dir="blastn_alignment"

# create output directory if it does not exist
mkdir -p "$output_dir"

# loop through each FASTA file in the input directory
for fasta_file in "$input_dir"/*.fas
do
  # get the base name of the file without the extension
  base_name=$(basename "$fasta_file" .fas)

  # define the output paths for the database and alignment
  db_output="$output_dir/${base_name}_db"
  alignment_output="$output_dir/${base_name}_alignment"

  # create the blast database 
  makeblastdb -in "$fasta_file" -out "$db_output" -parse_seqids -dbtype nucl

  # perform the blastn search
  blastn -query "$fasta_file" -db "$db_output" -word_size 4 -out "$alignment_output" -outfmt '6 qseqid sseqid pident'
done
