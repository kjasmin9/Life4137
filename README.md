# Mature miRNA analysis (sequence and co-occurrence in mammal species) 
by Jasmin Kaur
<br />
<br />
This repository contains all the code that has been used for the analysis of my Final Project
All code was developed by Jasmin Kaur unless stated otherwise.

## Part 1 - miRNA alignment and Co-occurrence distance analysis 
<br />
Refer to part_1/README_pt1.md for in detail tutorial of the analysis. 
<br />
part_1 directory contains all the scripts used in the analysis, and README_pt1.md provides an in-depth description of the order and how to run the scripts. <br />
<br />
The 'part_1/data files for part 1' directory contains FASTA files with their different versions obtained during the analysis for comparison purposes, in addition to a list of miRNAs to filter the FASTA files
<br />
<br />
<br />

## Part 2 - Ancestral reconstruction
<br />
The scripts for the analysis can be found in part2 directory. 
<br />
The directory contains three scripts, each script performs the same functions for the three different phenotypes: Placenta Shape, Membrane, and Venous Patter. 
<br />
<br />
The main tasks of each script are Ancestral reconstruction, miRNA analysis for divergent phenotypes, and statistical tests. 
<br />
The main input data are: <br />
1. A consensus tree in nexus format <br />
2. miRNA co-occurrence and phenotype data merged in one file (not available to the public) <br />
3. output from the part_1/dissimilarity_vs_co-occ.dist.R script: merged.df.pt1_1seq <br /> This file contains the dissimilarity values from aligned miRNA pairs.
<br />
<br />
<br />

##  Supplementary Data
<br />
Large files for supplementary material from the project have been included in the GitHub repository. Links for the files have been provided in the project write-up. 
