library(spaa) # packageVersion("spaa") ‘0.2.2’
getwd()


# with non-aligned miRNA pairs (manually inserted 1s) 
merged.df <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/1/merged.df_1seq", header = T) #from dissimilarity_vs_co-occ.dist.R, lines 134 & 135

#get different df for alignemnt scores and distance scores
al.dist <- merged.df[, c(2, 3, 5)]
mir.dist <- merged.df[, 2:4]

#to save
#write.csv(al.dist, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/1/al.dist_1seq", row.names = F,  col.names = T)
#write.csv(mir.dist, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/1/mir.dist_1seq",  row.names = F, col.names = T)

#run list2dist for each
al.dist1 <- list2dist(al.dist)
#list2dist for alignemnt is done!
mir.dist1 <- list2dist(mir.dist)
#list2dist for distance scores is done!

#convert both dists into matrices
al.dist.matrix <- as.matrix(al.dist1)
mir.dist.matrix <- as.matrix(mir.dist1)

#save the matrices!
write.csv(al.dist.matrix, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/1/al.dist.matrix_1seq", row.names = T,  col.names = T)
write.csv(mir.dist.matrix, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/1/mir.dist.matrix_1seq", row.names = T,  col.names = T)




#[2]
#with only aligned sequences
merged.df <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/na/merged.df.na_1seq", header = T)

#get different df for alignemnt scores and distance scores
al.dist <- merged.df[, c(2, 3, 5)]
mir.dist <- merged.df[, 2:4]

#to save
#write.csv(al.dist, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/na/al.dist_1seq", row.names = F,  col.names = T)
#write.csv(mir.dist, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/na/mir.dist_1seq",  row.names = F, col.names = T)

#run list2dist for each
al.dist1 <- list2dist(al.dist)
#list2dist for alignemnt is done!
mir.dist1 <- list2dist(mir.dist)
#list2dist for distance scores is done!

#convert both dists into matrices
al.dist.matrix <- as.matrix(al.dist1)
mir.dist.matrix <- as.matrix(mir.dist1)

#save the matrices!
write.csv(al.dist.matrix, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/na/al.dist.na.matrix_1seq", row.names = T,  col.names = T)
write.csv(mir.dist.matrix, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/na/mir.dist.na.matrix_1seq", row.names = T,  col.names = T) 





#by excluding na values dissimilarity (this reduces co-occurrence distance values) 
mirna_distance_dist.na <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/na/mir.dist.na.matrix_1seq", header = T)
alignment_dissimilarity_dist.na <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/na/al.dist.na.matrix_1seq", header = T)

rownames(mirna_distance_dist.na) <- mirna_distance_dist.na$X
mirna_distance_dist.na <- mirna_distance_dist.na[, 2:ncol(mirna_distance_dist.na)]

rownames(alignment_dissimilarity_dist.na) <- alignment_dissimilarity_dist.na$X
alignment_dissimilarity_dist.na <- alignment_dissimilarity_dist.na[, 2:ncol(alignment_dissimilarity_dist.na)]

alignment_dissimilarity_dist.na_omit <- as.dist(alignment_dissimilarity_dist.na)
alignment_dissimilarity_dist.na_omit <- dist2list(alignment_dissimilarity_dist.na_omit)

mirna_distance_dist.na_omit <- as.dist(mirna_distance_dist.na)
mirna_distance_dist.na_omit <- dist2list(mirna_distance_dist.na_omit)


mrg <- merge(alignment_dissimilarity_dist.na_omit, mirna_distance_dist.na_omit, by = c("col", "row"))
mrg <- na.omit(mrg)

alignment_dissimilarity_dist.na_omit<- mrg[, 1:3]
mirna_distance_dist.na_omit <- mrg[, c(1, 2, 4)]

al.dist.matrix <- list2dist(alignment_dissimilarity_dist.na_omit)
mir.dist.matrix <- list2dist(mirna_distance_dist.na_omit)

al.dist.matrix <- as.matrix(al.dist.matrix)
mir.dist.matrix <- as.matrix(mir.dist.matrix)

write.csv(al.dist.matrix, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/na.omit/al.dist.matrix_1seq", row.names = T,  col.names = T)
write.csv(mir.dist.matrix, file = "/gpfs01/home/mbxjk6/R/mantel/list2dist/na.omit/mir.dist.matrix_1seq", row.names = T,  col.names = T)
