# Part 1 analysis tutorial
 
This is a tutorial for miRNA analysis across mammal species.
Please follow the steps in chronological order for reproducibility purposes.

Main processes include: 
- Fasta file header string modification 
- Extraction of desired miRNA families
- Blastn alignment for multiple fasta files 
- Alignment similarity scores and miRNA co-occurrence analysis (Mantel test, Kruskal–Wallis test, Dunn's test & plotting) 


Versions of different fasta files and the list of miRNA families of interest are uploaded in Life4137/part_1/data files for part 1/ directory. The names of the files are kept the same as in the code.
The miRNA co-occurrence data file is part of an ongoing research project and is not available for public sharing.

All packages for RStudio can be installed by: install.packages('Package_name'), unless stated otherwise.


## get_fasta_files

In Bash: ```get_fasta_files```

This code contains the links for the 13 mamamals species FASTA files retrieved from MirGeneDB database. it concatenates all FASTA files into one file: mature_sequences.fas. The concatenated FASTA file can also be directly downloaded from: https://github.com/kjasmin9/Life4137/blob/main/part_1/data%20files%20for%20part%201/mature_sequences.fas



## get_miRNAs_of_interest.R

in RStudio: ```get_miRNAs_of_interest.R``` .
This code obtains a list of miRNAs of interest from the miRNA co-occurrence data frame
miRNAs are present as Column names. 

The list of miRNA families of interest (i.e., present in the co-occurrence distance matrix) can be obtained from: https://github.com/kjasmin9/Life4137/blob/main/part_1/data%20files%20for%20part%201/updated_mirnas_of.interest.csv



## Modify_fasta_headers.R

in RStudio: ```Modify_fasta_headers.R``` 
This code modifies the headers of mature_sequences.fas, by adding the the family name in the FASTA file identifiers. 
From ``` > Bta-Let-7-P1b_5p ```  to  ```>Let7 | Bta-Let-7-P1b_5p```.
Identifiers in this format are needed to filter the FASTA file for the miRNA families of interest using the seqkit command in bash. 




## Seqkit 
in bash: 
seqkit to get the miRNAs of interest.

To install Seqkit: https://github.com/shenwei356/seqkit/releases/tag/v2.8.2 
```
seqkit grep -f updated_mirnas_of.interest.csv modified_mature_sequences.fas > mirnas_of_interest.mature_sequences.fas
```


## create_single_csv.py

In bash: ```create_single_csv.py```
python script to obtain separate csv files, each containing a single miRNA family name (from the updated_mirnas_of.interest.csv). 
A part of this code was co-piloted with ChatGPT, an AI language model by OpenAI. 

**Usage**
```
chmod +x create_single_csv.py # give executable permissions
python3 create_single_csv.py # run the script 
```



## process_fasta.sh

In Bash: ```process_fasta.sh```
this script creates single fasta files for each miRNA family from the modified_mature_sequences.fas.
The script needs to be in the directory where csv outputs of create_single_csv.py are located

**Usage**
```
chmod +x process_fasta.sh
./process_fasta.sh
```



## Modify_fasta_headers_2.R

In RStudio: ```Modify_fasta_headers_2.R```
Modify headers of all .fas files present in the output directory from process_fasta.sh
This code allows for the fasta file identifiers to be compatible with makeblastdb. 
From ```>Let7 | Bta-Let-7-P1b_5p``` to ```> Bta-Let-7-P1b_5p | Let7```.
The input is a directory containing fasta files. 
This code was used to modify identifiers of the single miRNA family FASTA files and the FASTA file output from the seqkit command (mirnas_of_interest.mature_sequences.fas). The output of the second input FASTA file can be retrieved from: https://github.com/kjasmin9/Life4137/blob/main/part_1/data%20files%20for%20part%201/blast_modified_mature_sequences.fas

**Usage**: modify paths in the script



## blast.sh

In bash: ```blast.sh```
blastn: 2.14.1+ Package: blast 2.14.1, build Mar 13 2024 14:27:45
this script serves to run makeblastdb & blastn on all miRNA family .fas files.

Blast installation instructions: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html 


**Usage**
```
chmod +x blast.sh
sbatch blast.sh
```


## average_al.sh 

In bash: ```average_al.sh ```
This script calculates within miRNA family average similarity scores from blast alignment outputs.
A part of this code was co-piloted with ChatGPT, an AI language model by OpenAI 

**Usage**
```
chmod +x average_al.sh 
sbatch average_al.sh 
```


## extract_1st_seq.sh

In Bash: ```extract_1st_seq.sh```

This script extracts the first header and sequence from each fasta file in the input directory, the outputs are stored as a new fasta file
Since within-family similarity scores range from ~89% to 100%, one mature miRNA sequence is used to represent each miRNA family. 
The working directory for this code is **4blast_fas/** (created from) ```Modify_fasta_headers_2.R.```



## Blast alignment between miRNA families (all-vs-all)

In Bash: 
 
Run blastn for between miRNA families alignment (one sequence per family)

```
module load blast-uoneasy/2.14.1-gompi-2023a #activate blast if needed
makeblastdb -in all_fasta_1seq.fas -out db_output -parse_seqids -dbtype nucl
blastn -query all_fasta_1seq.fas -db db_output -word_size 4 -out alignment_1seq -outfmt '6 qseqid sseqid pident'
```

The code was adapted to run on the FASTA file containing all sequences for the 430 miRNA families:

```
module load blast-uoneasy/2.14.1-gompi-2023a #activate blast if needed
makeblastdb -in blast_modified_mature_sequences.fas -out db_output -parse_seqids -dbtype nucl
blastn -query blast_modified_mature_sequences.fas -db db_output -word_size 4 -out alignment4_mirnas_of_interest_wd4 -outfmt '6 qseqid sseqid pident'
````


## dissimilarity_vs_co-occ.dist.R

The ```dissimilarity_vs_co-occ.dist.R``` script contains all the analyses conducted on and between the alignment outputs from blastn searches and the miRNA co-occurrence data. 

Main analyses include: 

- miRNA string label modification
- matching miRNA pairs for alignment and mi-RNA co-occurrence data
- distance calculation for miRNA co-occurrence
- data frame merge based on ni-RNA pairs
- plots for distributions, regression, and boxplots
- Kruskal Wallis test, Dunn's test, and prop.test




## list2dist.R & sbatch_4_r.sh

This is an R script that converts a list object into dist object. 
The input file for this code is 'merged.df.pt1_1seq', obtained from the script ```dissimilarity_vs_co-occ.dist.R```.

```list2dist.R``` converts list objects of miRNA pairwise co-occurrence distance and dissimilarity scores into dist (matrix format) objects. 

```sbatch_4_r.sh``` serves to run ```list2dist.R``` as batch job in bash. 
list2dist (spaa package function), is time-consuming for large list objects, hence it is recommended to run the function as a background job.

**Usage** 
```
chmod +x list2dist.R
chmod +x sbatch_4_r.sh
sbatch sbatch_4_r.sh
```



## Run_mantel.R

The input for this script are the outputs from ```list2dist.R```
Mantel test tests for correlation between two matrices. Pairwise miRNA co-occurrence distance and pairwise miRNA alignment dissimilarity values stored each as dist objects can be correlated. 

**Usage**: run ```Run_mantel.R``` in RStudio.










## References
- Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.
- Wei Shen, Shuai Le, Yan Li*, and Fuquan Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.
- OpenAI, 2024. *ChatGPT: OpenAI's GPT-4 Language Model*. Available at: https://www.openai.com/chatgpt [Accessed 2 July 2024].

