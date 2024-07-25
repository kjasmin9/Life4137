#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50g
#SBATCH --time=24:00:00
#SBATCH --job-name=averaging


# dfine input and output directories
input_dir="Path/to/blast.sh/output_dir/all_alignments"
output_dir="al_averages"

# create output directory if it does not exist
mkdir -p "$output_dir"

# loop through each alignment file in the input directory
for alignment in "$input_dir"/*_alignment
do
  # get the base name of the file without the extension
  base_name=$(basename "$alignment" _alignment)

  # define the output path for averages of each alignment file
  average_output="$output_dir/${base_name}_avg_al"

  # average the alignments for each file and save the result to the output file (with miRNA family name)
  # This code was co-piloted with ChatGPT, an AI language model by OpenAI.
  awk '{ sum += $3; count++ } END { if (count > 0) print sum / count; }' "$alignment" > "$average_output"

done
