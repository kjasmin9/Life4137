'''if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")'''
#BiocManager::install("Biostrings")
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
