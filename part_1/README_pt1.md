# Part 1 analysis tutorial
 
This is a tutorial for miRNA analysis across mammal species.
Please follow the steps in chronological order for reproducibility purposes.

Main processes include: 
- Fasta file header string modification 
- Extraction of desired miRNA families
- Blastn alignment for multiple fasta files 
- Alignment similarity scores and miRNA co-occurrence analysis (Mantel test, Kruskal–Wallis test, Dunn's test & plotting) 


Versions of different fasta files and the list of miRNA families of interest are uploaded in Life4137/part_1/data files for part 1/ directory. The names of the files are kept the same as in the code.
The miRNA co-occurrence data file is part of an ongoing research project and not available for public sharing.

All packages for RStudio can be installed by: install.packages('Package_name'), unless stated otherwise.


## get_fasta_files

In Bash: ```get_fasta_files```

```
# download the first 13 mammal mature miRNA seqs from https://mirgenedb.org/download exclude the all_species (since it has other than mammals) 
mkdir fasta_files 
wget https://mirgenedb.org/fasta/hsa?mat=1
wget https://mirgenedb.org/fasta/mml?mat=1
wget https://mirgenedb.org/fasta/mmu?mat=1
wget https://mirgenedb.org/fasta/rno?mat=1
wget https://mirgenedb.org/fasta/cpo?mat=1
wget https://mirgenedb.org/fasta/ocu?mat=1
wget https://mirgenedb.org/fasta/cfa?mat=1
wget https://mirgenedb.org/fasta/bta?mat=1
wget https://mirgenedb.org/fasta/dno?mat=1
wget https://mirgenedb.org/fasta/ete?mat=1
wget https://mirgenedb.org/fasta/sha?mat=1
wget https://mirgenedb.org/fasta/mdo?mat=1
wget https://mirgenedb.org/fasta/oan?mat=1

# merge the fasta files into one
cat *.fas > mature_sequences.fas
```

## get_miRNAs_of_interest.R

in RStudio: ```get_miRNAs_of_interest.R``` .
This code obtains a list of miRNAs of interest from the miRNA co-occurrence data frame
miRNAs are present as Column names

```
#read in the data
mirnas_of_interest <- read.csv('/Path/to/updated_df1.csv', header = F )
mirnas_of_interest <- mirnas_of_interest[1, c(5:ncol(mirnas_of_interest))] # remove phenotype columns
mirnas_of_interest <- t(mirnas_of_interest)
#insert '-' afer miRNA family name prefix 
mirnas_of_interest1 <- gsub("^Mir", "Mir-", mirnas_of_interest[, 1])
mirnas_of_interest1 <- gsub("^Novel", "Novel-", mirnas_of_interest1)
mirnas_of_interest1 <- gsub("^Let", "Let-", mirnas_of_interest1)
mirnas_of_interest1 <- gsub("^Bantam", "Bantam-", mirnas_of_interest1)
mirnas_of_interest1 <- as.matrix(mirnas_of_interest1 )
#saving mirnas_of_interest to filter out the fasta file 
write.csv(mirnas_of_interest1, file = "/Desired/path/to/mirnas_of_interest1.csv", row.names = FALSE, col.names = FALSE)
```

## Modify_fasta_headers.R

in RStudio: ```Modify_fasta_headers.R``` 
This code modifies the headers of mature_sequences.fas

```
'''if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")'''
#BiocManager::install("Biostrings")
library('Biostrings')
packageVersion("Biostrings") # ‘2.70.3’

fasta <- readDNAStringSet("Path/to/mature_sequences.fas")

# Get original and modified names
og_names <- names(fasta)

#define functions for modifications
modify_name1 <- function(name) {
  sub(".*?-", "", name)
}

modify_name2 <- function(name) {
  sub("^(([^-]*-){1}[^-]*)-.*", "\\1", name)
}

modify_name3 <- function(name) {
  sub("_.*", "", name)
}

modify_name4 <- function(name) {
  gsub("-", "", name)
}

# apply name modification functions
modified_names <- sapply(og_names, modify_name1)
modified_names <- sapply(modified_names, modify_name2)
modified_names <- sapply(modified_names, modify_name3)
modified_names <- sapply(modified_names, modify_name4)

# Combine original and modified names, separate by pipe '|' or desired separator 
new_names <- paste(modified_names, og_names, sep=" | ")

#assign the updated names to the fasta file 
names(fasta) <- new_names
# Write modified Fasta file
writeXStringSet(fasta, "PAth/to/modified_mature_sequences.fas", format = "fasta")
```


## Seqkit 
in bash: 
seqkit to get the miRNAs of interest 
to install Seqkit: https://github.com/shenwei356/seqkit/releases/tag/v2.8.2 
```
seqkit grep -f updated_mirnas_of.interest.csv modified_mature_sequences.fas > mirnas_of_interest.mature_sequences.fas
```


## create_single_csv.py

In bash: ```create_single_csv.py```
python script to obtain separate csv files, each containinag a single miRNA family name (from the updated_mirnas_of.interest.csv)

```
#!/bin/bash

import csv
import os

# define the path to the input CSV file
input_csv_path = 'Path/to/updated_mirnas_of.interest.csv'

# define and create the directory to store the output files
output_dir = 'Desired/path/to/mirnas_of_interest_482''
os.makedirs(output_dir, exist_ok=True)

# read the input CSV file
with open(input_csv_path, newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        content = row[0]
        filename = f"{content}.csv"
        file_path = os.path.join(output_dir, filename)

        # write the content to the new file
        with open(file_path, 'w', newline='') as newfile:
            writer = csv.writer(newfile)
            writer.writerow([content])

print(f'Written .csv file can be found in {output_dir}')
```
Usage
```
chmod +x create_single_csv.py # give executable permissions
python3 create_single_csv.py # run the script 
```

## process_fasta.sh

In Bash: ```process_fasta.sh```
this is a script to create single fasta files for each miRNA family from the modified_mature_sequences.fas
The script needs to be in the directory where csv outputs of create_single_csv.py are located
nano process_fasta.sh

```
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
```
Usage: 
```
chmod +x process_fasta.sh
./process_fasta.sh
```


## Modify_fasta_headers_2.R

In RStudio: ```Modify_fasta_headers_2.R```
Modify headers of all .fas files present in the output directory from process_fasta.sh
This code allows for the fasta file identifiers to be compatible with makeblastdb
the input is a directory containing fasta files. 
Usage: modify paths in the script

```
library(Biostrings) # packageVersion("Biostrings") ‘2.70.3’
library(foreach) # packageVersion('foreach') ‘1.5.2’
library(doParallel) # packageVersion('doParallel') ‘1.0.17’

# To run functions in a parallel way --> code runs faster
numCores <- detectCores() - 1  # Use one less than the total number of cores available
cl <- makeCluster(numCores)
registerDoParallel(cl)


#set input and output directories (create output dir if it deos not exist)
input_dir <- "C:/Path/to/only_fas/"
output_dir <- "C:/PAth/to/4blast_fas/"
dir.create(output_dir)

#store all the fasta files paths in a list
fasta_files <- list.files(input_dir, pattern = "\\.fas", full.names = T)

#now run code on all files at once
foreach(fasta_file = fasta_files, .packages = "Biostrings") %dopar% {
  fasta <- readDNAStringSet(fasta_file)
  
  fasta_ids <- names(fasta) #store fasta ids 
  #apply desired cahnges to the fasta ids
  fam_name <- sub("^([^|]+)\\|.+$", "\\1", fasta_ids)
  member_name <- sub("^.+\\|(.*)$", "\\1", fasta_ids)
  #merge the names in in the identifiers so that miRNA name precedes miRNA family prefix
  adjusted_names <- paste(member_name, fam_name, sep = " | ")
  
  #assign new fasta ids to the fasta file
  names(fasta) <- adjusted_names
  #assign nemas to the output files 
  output_file <- file.path(output_dir, basename(fasta_file))
  
  writeXStringSet(fasta, output_file, format = "fasta")
}

stopCluster(cl)
```


## blast.sh

In bash: ```blast.sh```
blastn: 2.14.1+ Package: blast 2.14.1, build Mar 13 2024 14:27:45
this script serves to run makeblastdb & blastn on all miRNA family .fas files
Blast installation instructions: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html 

```
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
```
Usage
```
chmod +x blast.sh
sbatch blast.sh
```


## average_al.sh 

In bash: ```average_al.sh ```
This script calculates within miRNA family average similarity scores from blast alignemnt outputs.

```
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
  awk '{ sum += $3; count++ } END { if (count > 0) print sum / count; }' "$alignment" > "$average_output"

done
```
Usage
```
chmod +x average_al.sh 
sbatch average_al.sh 
```


## extract_1st_seq.sh

In Bash: ```extract_1st_seq.sh```

This script extracts the first header and sequence from each fasta file in the input directory, the outputs are stored as a new fasta file
Since within-family similarity scores range from ~89% to 100%, one mature miRNA sequence is used to represent each miRNA family. 
The working directory for this code is **4blast_fas/** (created from) Modify_fasta_headers_2.R.

```
#!/bin/bash

# create empty output file
output_file="all_fasta_1seq.fas" > $output_file

# loop through each fasta file in the directory and extract first ID and sequence
for file in *.fas
do
  # output in output_file
  head -n 2 "$file" >> $output_file
done
```
Usage
```
chmod +x extract_1st_seq.sh
./extract_1st_seq.sh 
```

## Blast alignment between miRNA families (all-vs-all)

In Bash: 
 
Run blastn for between miRNA families alignment

```
module load blast-uoneasy/2.14.1-gompi-2023a
makeblastdb -in all_fasta_1seq.fas -out db_output -parse_seqids -dbtype nucl
blastn -query all_fasta_1seq.fas -db db_output -word_size 4 -out alignment_1seq -outfmt '6 qseqid sseqid pident'
```



insert the dissimilarity_vs_co-occ.dist.R here



## list2dist.R & sbatch_4_r.sh
In bash write an R script to convert a list object into dist object.
Use the scripts ```list2dist.R``` and ```sbatch_4_r.sh```.
```list2dist.R``` converts list objects of miRNA pairwise co-occurrence distance and dissimilarity scores into dist (matrix format) objects. 
```sbatch_4_r.sh``` serves to run ```list2dist.R``` as batch job. 
list2dist (spaa package function), is time-consuming for large list objects, hence it is recommended to run the function as a background job.

Usage 
```
chmod +x list2dist.R
chmod +x sbatch_4_r.sh
sbatch sbatch_4_r.sh
```






## References
- Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.
- Wei Shen, Shuai Le, Yan Li*, and Fuquan Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.
