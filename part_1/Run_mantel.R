library(spaa) # packageVersion("spaa") ‘0.2.2’
library(vegan) #  packageVersion("vegan") ‘2.6.2’

######~na.omit~######
mirna_distance_dist.na.omit <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/na.omit/mir.dist.matrix_1seq", header = T)
alignment_dissimilarity_dist.na.omit <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/na.omit/al.dist.matrix_1seq", header = T)

# get dist data ready for mantel test
# assign miRNA names as rownames for the data.frames
rownames(mirna_distance_dist.na.omit) <- mirna_distance_dist.na.omit$X
mirna_distance_dist.na.omit <- mirna_distance_dist.na.omit[, 2:ncol(mirna_distance_dist.na.omit)]

# assign miRNA names as row names
rownames(alignment_dissimilarity_dist.na.omit) <- alignment_dissimilarity_dist.na.omit$X
alignment_dissimilarity_dist.na.omit <- alignment_dissimilarity_dist.na.omit[, 2:ncol(alignment_dissimilarity_dist.na.omit)]

# Perform mantel test
mantel.test.na.omit <- mantel(mirna_distance_dist.na.omit, alignment_dissimilarity_dist.na.omit, method = "spearman", permutations = 9999, na.rm = T)
print(mantel.test.na.omit)



######~na~######
mirna_distance_dist.na <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/na/mir.dist.na.matrix_1seq", header = T)
alignment_dissimilarity_dist.na <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/na/al.dist.na.matrix_1seq", header = T)

#get dist data ready for mantel test
#assign miRNA names as rownames for the data.frames
rownames(mirna_distance_dist.na) <- mirna_distance_dist.na$X
mirna_distance_dist.na <- mirna_distance_dist.na[, 2:ncol(mirna_distance_dist.na)]

# assign miRNA names as row names
rownames(alignment_dissimilarity_dist.na) <- alignment_dissimilarity_dist.na$X
alignment_dissimilarity_dist.na <- alignment_dissimilarity_dist.na[, 2:ncol(alignment_dissimilarity_dist.na)]

# Perform mantel test
mantel.test.na <- mantel(mirna_distance_dist.na, alignment_dissimilarity_dist.na, method = "spearman", permutations = 9999, na.rm = T)
print(mantel.test.na)



######~1~######
mirna_distance_dist.1 <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/1/mir.dist.matrix_1seq", header = T)
alignment_dissimilarity_dist.1 <- read.csv("/gpfs01/home/mbxjk6/R/mantel/list2dist/1/al.dist.matrix_1seq", header = T)

#get dist data ready for mantel test
#assign miRNA names as rownames for the data.frames
rownames(mirna_distance_dist.1) <- mirna_distance_dist.1$X
mirna_distance_dist.1 <- mirna_distance_dist.1[, 2:ncol(mirna_distance_dist.1)]

# assign miRNA names as row names
rownames(alignment_dissimilarity_dist.1) <- alignment_dissimilarity_dist.1$X
alignment_dissimilarity_dist.1 <- alignment_dissimilarity_dist.1[, 2:ncol(alignment_dissimilarity_dist.1)]

# Perform mantel test
mantel.test.1 <- mantel(mirna_distance_dist.1, alignment_dissimilarity_dist.1, method = "spearman", permutations = 9999, na.rm = T)
print(mantel.test.1)
