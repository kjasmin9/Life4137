library(phytools) # packageVersion('phytools') ‘2.1.1’
library(readxl) # packageVersion('readxl')  ‘1.4.3’
library(ape) # packageVersion('ape') ‘5.8’
library(tidyverse) # packageVersion('tidyverse') ‘2.0.0’
library(FSA) #packageVersion('FSA') ‘0.9.5’
library(RColorBrewer) # packageVersion('RColorBrewer') ‘1.1.3
library(ggplot2) # packageVersion('ggplot2') ‘3.4.4’



# read in consensus tree & phenotype data
tree <- read.tree('C:/Users/jaxk2/Downloads/consensus_tree.phy.txt') 
df <- read_xlsx('C:/Users/jaxk2/Downloads/updated_df.xlsx', sheet =  'Sheet2')


#[A]  ancestral reconstruction 
##############################

any(df$`Placental Shape` == 'NA') #TRUE,  therefore remove the NAs
df.placenta <- df[df$`Placental Shape` != "NA", ] #remove the NAs from the df$`Placenta Shape`

# get common species between tree and df 
placenta.shared_species <- intersect(df.placenta$Species, tree$tip.label)

#retain only common species in both tree and df
tree.placenta <- keep.tip(tree, placenta.shared_species)
df.placenta <- df.placenta[df.placenta$Species %in% placenta.shared_species,]

#reorder species in the df according to tree
try <- match(tree.placenta$tip.label, df.placenta$Species)
df.placenta <- df.placenta[try, ]
df.placenta <- df.placenta[, c(1, 4)] #df with only species and Placental_shape columns

#add phenotype onto tree as states
tree.placenta$states <- as.character(df.placenta$`Placental Shape`)



# following the tutorial http://www.phytools.org/eqg2015/asr.html
## looking at the whole distribution from a sample of stochastic maps
### ancestral reconstruction with creation of stochastic trees
df.placenta1 <- as.data.frame(df.placenta)                                      
rownames(df.placenta1) <- df.placenta$Species

fmode <- setNames(df.placenta1$`Placental Shape`, rownames(df.placenta1))
fmode1 <- as.factor(fmode)

FMODE<-to.matrix(fmode1,levels(fmode1))
FMODE<-FMODE[tree.placenta$tip.label,] #to make sure that species in FMODE2 and the tree match 
x <- fmode1 #easier for code writing

dff <- brewer.pal(4, "Set1") #color-Blind friendly palette (4 for num of phenotypes)
#create color palette based on number of phenotypes
cols <- setNames(palette()[1:length(unique(x))],sort(unique(x))) #remove dff if you want to use palette only

# generate 100 stochastic character trees (model = ER assumes equal probs of each character being ancestral)
mtrees<-make.simmap(tree.placenta,x,model="ER",nsim=100)
pd<-summary(mtrees,plot=FALSE) #summary of the 100 stochastic trees
plot(pd,fsize=0.6,ftype="i") #to visualise the summary of the stochastic trees

plot(mtrees[[1]],cols,fsize=0.7,ftype="i")
tiplabels(pie=FMODE,piecol=palette(cols),cex=0.3)
nodelabels(pie=pd$ace,piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=F,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree.placenta)),fsize=1)




#[B] getting unique miRNAs for divergently evolved phenotype groups
###################################################################
group1 <- 'Otolemur crassicaudatus 
Daubentonia madagascariensis 
Indri indri 
Microcebus murinus 
Varecia variegata 
Lemur catta
Eulemur fulvus'

# from list to data.frame
group1 <- as.data.frame(strsplit(group1, "\n"))
group1$Species <- sub(" ", "_", group1$c..Otolemur.crassicaudatus.....Daubentonia.madagascariensis....) #,match species names to names of species in df
group1$Species <- sub(" ", "", group1$Species)
group1 <- group1 %>% select(Species)
#retain only species present in phenotypes data frame
species1 <- df %>% filter(df$Species %in% group1$Species)
species1 <- species1[, c(1, 5:ncol(species1))]

#calculate mean for each mirna family (percentage occurrence)
means <- colMeans(species1[, 2:ncol(species1)], na.rm = T) # calculate means, results will be in new df)
means_row <- data.frame(occurrence.perc = means) #
means_row <- means_row %>% add_row(.before = 1) #adding a row to match the num of cols in species1 
rownames(means_row) <- colnames(species1) #match row and col names of the two dfs before merging 
means_row[1, 1]<- "occ_perc"
species1 <- rbind(species1, t(means_row)) # add the new means_row into the species1 df

#change species col into rownames
species1 <- as.data.frame(species1)
rownames(species1) <- species1$Species
species1 <- species1[, c(2:ncol(species1))]

##################################################################################MAIN STEP
#filter only mirnas (cols), which have occurrence == 100% across all species
species1.100mirna <- as.data.frame(t(species1))
species1.100mirna <- species1.100mirna %>% filter(occ_perc == 1)
species1_list100 <- as.list(rownames(species1.100mirna)) # save names of mirnas occurring in species 1 (with occurrence100%), as a list
# ↑↑↑↑↑ use space above to get different occurrence percentages i.e, >=90% ↑↑↑↑↑
#################################################################################MAIN STEP

# repeat all that for group2 species 

group2 =  "Catagonus wagneri
Sus scrofa
Tragulus javanicus
Hippopotamus amphibius
Megaptera novaeangliae
Cephalorhynchus commersonii
Orcinus orca
Phocoena phocoena
Delphinapterus leucas
Inia geoffrensis
Camelus bactrianus
Lama glama
Diceros bicornis
Ceratotherium simum
Rhinoceros unicornis
Equus asinus"

group2 <- as.data.frame(strsplit(group2, "\n"))
group2$Species <- sub(" ", "_", group2$c..Catagonus.wagneri....Sus.scrofa....Tragulus.javanicus....Hippopotamus.amphibius...)
group2$Species <- sub(" ", "", group2$Species)
group2 <- group2 %>% select(Species)
species2 <- df %>% filter(df$Species %in% group2$Species)
species2 <- species2[, c(1, 5:ncol(species2))]


means2 <- colMeans(species2[, 2:ncol(species2)], na.rm = T) # calculate means, results will be in new df)
means_row2 <- data.frame(occurrence.perc = means2) #
means_row2 <- means_row2 %>% add_row(.before = 1) #adding a row to match the num of cols in species1 
rownames(means_row2) <- colnames(species2) #match row and col names of the two dfs before merging 
means_row2[1, 1]<- "occ_perc"
species2 <- rbind(species2, t(means_row2)) # add the new means_row into the species1 df

#change species col into rownames
species2 <- as.data.frame(species2)
rownames(species2) <- species2$Species
species2 <- species2[, c(2:ncol(species2))]

#filter only mirnas (cols), which have occurrence == 100% across all species
species2.100mirna <- as.data.frame(t(species2))
species2.100mirna <- species2.100mirna %>% filter(occ_perc == 1)
species2_list100 <- as.list(rownames(species2.100mirna)) # save names of mirnas occurring in species 1 (with occurrence100%), as a list


#[2]
# obtain: mirnas in 1 but not in 2, and mirnas in 2, but not in 1
## unique.species1.100occ <- mirnas unique to group1, with 100% occurrence across all species in group 1
## unique.species2.100occ <- mirnas unique to group2, with 100% occurrence across all species in group 2
## similar.species1.2.100occ <- both mirnas present in group 1 and group 2

unique.species1.100occ <- setdiff(species1_list100, species2_list100)
unique.species2.100occ <- setdiff(species2_list100, species1_list100)
similar.species1.2.100occ <- intersect(species2_list100, species1_list100)
'similar.species1.2.100occ1 <- intersect(species1_list100, species2_list100) 
remove(similar.species1.2.100occ1)' #they are the same
#Mir1247 is in unique.species1.100occ
any(rownames(species1.100mirna) == 'Mir1247') #checking #TRUE
any(rownames(species2.100mirna) == 'Mir1247') #checking #FALSE
# checked: this means you named the unique mirnas correctly, ie unique.species1.100occ belong to species in group 1. 

#save the unique lists
#some modifications before saving 
unique.species1.100occ <- t(as.data.frame(unique.species1.100occ))
unique.species2.100occ <- t(as.data.frame(unique.species2.100occ))
similar.species1.2.100occ <- t(as.data.frame(similar.species1.2.100occ))

#save!  
'write.csv(unique.species1.100occ, file = "C:/Users/jaxk2/Downloads/unique.species1.100occ", row.names = F , quote = F)
write.csv(unique.species2.100occ, file = "C:/Users/jaxk2/Downloads/unique.species2.100occ", row.names = F, quote = F )
write.csv(similar.species1.2.100occ, file = "C:/Users/jaxk2/Downloads/similar.species1.2.100occ", row.names = F, quote = F)'

#end result: three lists of miRNAs (unique to gorup1, unique to group2, shared between the two gorups)



#[C] there are miRNA lists for each group, so how do the alignments between these groups look like?
# this code combines each group mirnas with each other, creating all possible pairings between the two 
# for each pair, mature miRNA alignement values (dissimilarity from part 1) are obtained
# 3 alignment values are then compared with Krusk-Wallis test, and Dunn test (if Krusk Wallis is significant)
#################################################################################################
grp1 <- read.csv("C:/Users/jaxk2/Downloads/unique.species1.100occ")
grp2 <- read.csv("C:/Users/jaxk2/Downloads/unique.species2.100occ")
shared <-read.csv("C:/Users/jaxk2/Downloads/similar.species1.2.100occ")
#read in the pairwise mirna co-occurrence distance and dissimilarity average lists 
merged.df.pt1 <- read.csv("C:/Users/jaxk2/Downloads/merged.df.pt1_1seq", header = T)
merged.df.pt1 <- merged.df.pt1 %>% select(merged_mirna, distance_score, dissimilarity_average, mirna1.x, mirna2.x)
merged.df.pt1 <- merged.df.pt1[merged.df.pt1$dissimilarity_average != 1, ] #dsicard manually inserted 1s (for non-aligned pairs)


#[1] mutually exclusive mirnas group1 vs group2 and dissimilarity and distance scores. (although main analysis is on dissimilarity scores)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add grp2 as the top row
g1.vs.g2 <- grp1 %>% add_column(t(grp2)) #this adds grp2 as the first row (each cell value is repeated in the column)
grp2_ <- as_data_frame(g1.vs.g2$`t(grp2)`) # converts the added row into data.frame (previously matrix)
g1.vs.g2 <- bind_cols(grp1$V1, grp2_) # now dim: 40x6 (prev it was 40x2)

# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'a <- g1.vs.g2 %>% select(...1, V1)
a1 <- g1.vs.g2 %>% select(...1, V2)
a2 <- g1.vs.g2 %>% select(...1, V3)
a3 <- g1.vs.g2 %>% select(...1, V4)
a4 <- g1.vs.g2 %>% select(...1, V5)
df.list <- list(a, a1, a2, a3, a4)'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df.list1 <- lapply(seq(2, ncol(g1.vs.g2)), function(i) {
  bind_cols(grp1$V1, g1.vs.g2[, i])
}) #this command does all the above in the string, without being repetitive 

#change colnames for each df
df.list1 <- lapply(df.list1, function(x) {
  colnames(x) <- c("group1", "group2")
  return(x)
})


#merge vertically each df (in df.list) into one
g1.vs.g2 <- reduce(df.list1, rbind) # same code as: rbind(df.list[[1]], df.list[[2]], df.list[[3]], df.list[[4]], df.list[[5]])
g1.vs.g2 <- g1.vs.g2 %>% unite(merged_mirna, group1, group2, sep = '-')# merge the mirna names to match the merged.df.pt1


#merge g1.vs.g2 and merged.df.pt1 
g1.vs.g2_merged.df.pt1 <- merge(g1.vs.g2, merged.df.pt1, by = "merged_mirna", all.x = T)
g1.vs.g2_merged.df.pt1 <- na.omit(g1.vs.g2_merged.df.pt1) #any(merged.df.pt1 == "Mir3613") # FALSE, tat means that the mir3613 was not present in updated.df 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#[2] group1 mutually exclusive vs shared mirnas
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g1.vs.shared <- grp1 %>% add_column(t(shared))
shared_ <- as_data_frame(g1.vs.shared$`t(shared)`)
g1.vs.shared <- bind_cols(grp1$V1, shared_)


df.list_1sh <- lapply(seq(2, ncol(g1.vs.shared)), function(i) {
  bind_cols(grp1$V1, g1.vs.shared[, i])
})

df.list_1sh <- lapply(df.list_1sh, function(x) {
  colnames(x) <- c("group1", "shared_group")
  return(x)
})

#merge vertically each df in the list
g1.vs.shared <- reduce(df.list_1sh , rbind)
g1.vs.shared <- g1.vs.shared %>% unite(merged_mirna, group1, shared_group, sep = '-')# merge the mirna names to match the merged.df.pt1
#merge g1.vs.g2 and merged.df.pt1 
g1.vs.shared_merged.df.pt1 <- merge(g1.vs.shared, merged.df.pt1, by = "merged_mirna", all.x = T)
g1.vs.shared_merged.df.pt1 <- na.omit(g1.vs.shared_merged.df.pt1) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#[3] group2 mutually exclusive vs shared mirnas
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

g2.vs.shared <- grp2 %>% add_column(t(shared))
shared_2 <- as_data_frame(g2.vs.shared$`t(shared)`)
g2.vs.shared <- bind_cols(grp2$V1, shared_2)


df.list_2sh <- lapply(seq(2, ncol(g2.vs.shared)), function(i) {
  bind_cols(grp2$V1, g2.vs.shared[, i])
})

df.list_2sh <- lapply(df.list_2sh, function(x) {
  colnames(x) <- c("group2", "shared_group")
  return(x)
})

#merge vertically each df in the list
g2.vs.shared <- reduce(df.list_2sh , rbind)
g2.vs.shared <- g2.vs.shared %>% unite(merged_mirna, group2, shared_group, sep = '-')# merge the mirna names to match the merged.df.pt1


g2.vs.shared_merged.df.pt1 <- merge(g2.vs.shared, merged.df.pt1, by = "merged_mirna", all.x = T)
g2.vs.shared_merged.df.pt1 <- na.omit(g2.vs.shared_merged.df.pt1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#[4] plot the dissimilarity values for the 3 different groups
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#add label cols to all the 3 dfs ##merge vertically the 3 dfs ($dissimilarity & $label cols) 
#g1.vs.g2_merged.df.pt1 #g1.vs.shared_merged.df.pt1 #g2.vs.shared_merged.df.pt1

g1.vs.g2_merged.df.pt1 <- g1.vs.g2_merged.df.pt1 %>% 
  add_column(label = rep("g1.vs.g2", nrow(g1.vs.g2_merged.df.pt1)))

g1.vs.shared_merged.df.pt1 <- g1.vs.shared_merged.df.pt1 %>% 
  add_column(label = rep("g1.vs.shared", nrow(g1.vs.shared_merged.df.pt1)))

g2.vs.shared_merged.df.pt1 <- g2.vs.shared_merged.df.pt1 %>% 
  add_column(label = rep("g2.vs.shared", nrow(g2.vs.shared_merged.df.pt1)))
#merge vertically dissimilarity values of each df along with the labels (categorical groups) 
dissimilarity_3groups <- rbind(g1.vs.g2_merged.df.pt1[, c(3, 6)],g1.vs.shared_merged.df.pt1[, c(3, 6)],g2.vs.shared_merged.df.pt1[, c(3, 6)] )


#~~~~~~~~~~~~~~~~~~~~~ for KW & DUNN
#BOXPLOTS
ggplot(dissimilarity_3groups) + 
  aes(x = label, y = dissimilarity_average, fill = label) + 
  geom_boxplot()

#see dissimilarity distribution for each group
##first group_by the categorical variables (i.e., the 3 groups) 
dissimilarity_3groups <- dissimilarity_3groups %>%
  group_by(label)

ll <- dissimilarity_3groups %>%
  ggplot(aes(x = dissimilarity_average)) +
  geom_histogram(aes(fill = label), bins = 25) +
  facet_wrap(~ label) +
  labs(title = "miRNA Dissimilarity distribution (Placental Shape)", 
       x = "miRNA sequence Dissimilarity", 
       fill = 'Groups') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) #plot title in the center
  

#ggsave("C:/Users/jaxk2/Downloads/Placenta_dissimilarity_distribution_supp.png", plot = ll, width = 6, height = 4, units = "in", dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~ for proportion of zeros in each group  
#proportion of 0s in each group
dissimilarity_3groups <- dissimilarity_3groups <- dissimilarity_3groups %>%
  mutate(is_zero = dissimilarity_average == 0)

proportions <-  dissimilarity_3groups %>%
  group_by(label, is_zero) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

#plot the proportions
kk <- ggplot(proportions, aes(x = label, y = proportion, fill = as.factor(is_zero))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("blue", "orange"), labels = c("Identical", "Non-Identical")) +
  labs(x = "Comparison Groups", y = "Proportion", fill = " ", 
       title = 'Proportion of Identical and non-Identical Sequences') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("C:/Users/jaxk2/Downloads/Placenta_proportions_MAIN.png", plot = kk, width = 6, height = 4, units = "in", dpi = 300)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#[5] Kruskal and DUNN's test 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
krusk <- kruskal.test(dissimilarity_average ~ label, data = dissimilarity_3groups)
dunn <- dunnTest(dissimilarity_average ~ label, data = dissimilarity_3groups, method = "bonferroni")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#[6] proportions test
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g1g2_vs_g1shared_test <- prop.test(c(20, 568), c(9+20, 160+568)) #X-squared = 0.84854, df = 1, p-value = 0.357
g1g2_vs_g2shared_test <- prop.test(c(20, 85), c(9+20, 19+85)) #X-squared = 1.5216, df = 1, p-value = 0.2174
g1shared_vs_g2shared_test <- prop.test(c(568, 85), c(160+568, 19+85)) #X-squared = 0.53792, df = 1, p-value = 0.4633

#none is significant
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

