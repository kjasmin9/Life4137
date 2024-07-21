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
