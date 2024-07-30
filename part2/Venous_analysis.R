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

any(df$`Venous Pattern` == 'NA') #TRUE thf remove the NAs
df.venous <- df[df$`Venous Pattern` != "NA", ] #remove the NAs from the df$`Venous Pattern` #only 97 out of 300 species have membrance phenotype data
any(df.venous$`Venous Pattern` == 'NA') #check, it is to be FALSE

# get common species between tree and df 
venous.shared_species <- intersect(df.venous$Species, tree$tip.label)

#retain only common species in both tree and df
tree.venous <- keep.tip(tree, venous.shared_species)
df.venous <- df.venous[df.venous$Species %in% venous.shared_species,]

#reorder species in the df according to tree
try <- match(tree.venous$tip.label, df.venous$Species)
df.venous <- df.venous[try, ]
df.venous <- df.venous[, c(1, 3)] #df with only species and VEnous patterns columns

#add phenotype onto tree as states
tree.venous$states <- as.character(df.venous$`Venous Pattern`)



# following the tutorial http://www.phytools.org/eqg2015/asr.html
## looking at the whole distribution from a sample of stochastic maps
### ancestral reconstruction with creation of stochastic trees
df.venous1 <- as.data.frame(df.venous)
rownames(df.venous1) <- df.venous1$Species

fmode3 <- setNames(df.venous1$`Venous Pattern`, rownames(df.venous1))
fmode3 <- as.factor(fmode3)

#need this for make.simmap
FMODE3<-to.matrix(fmode3,levels(fmode3))
FMODE3<-FMODE3[tree.venous$tip.label,] #to make sure that species in FMODE3 and the tree match 

# get list of phenotyoes
unique(df.venous1$`Venous Pattern`) #[1] "Labyrinthine" "Villous" "Trabecular" : 3 phenotypes, hence 3 cols for the palette
dff3 <- brewer.pal(3, "Set1") #color-Blind friendly palette (3 for num of phenotypes)
#create color palette based on number of phenotypes
cols <- setNames(palette(dff3)[1:length(unique(fmode3))],sort(unique(fmode3)))




# generate 100 stochastic character trees (model = ER assumes equal probs of each character being ancestral)
mtrees.venous <- make.simmap(tree.venous, fmode3, model = "ER", nsim = 100)
pd.venous <-summary(mtrees.venous,plot=FALSE) #summary of the 100 stochastic trees
plot(pd.venous,fsize=0.5,ftype="i")#to visualise the summary of the stochastic trees

#plot a stochastic tree with ancestral state probabilities (states at the nodes)
plot(mtrees.venous[[1]],cols,fsize=0.7,ftype="i")
tiplabels(pie=FMODE3,piecol=palette(cols),cex=0.34) 
nodelabels(pie=pd.venous$ace,piecol=cols,cex=0.34) 
add.simmap.legend(colors=cols,prompt=T,x=0.5*par()$usr[1],
                  y=-max(nodeHeights(tree.venous)),fsize=0.99)



#[B] getting unique miRNAs for divergently evolved phenotype groups
###################################################################
group1.v <- 'Colobus angolensis
Pygathrix nemaeus
Nasalis larvatus
Semnopithecus entellus
Vacaca mulatta
Mandrillus sphinx
Miopithecus talapoin
Erythrocebus patas
an troglodytes
Pan paniscus lomo sapiens
Gorilla gorilla
Eulemur fulvus
Lemur catta
Varecia variegata
Microcebus murinus
Indri indri
Daubentonia madagascariensis
Otolemur crassicaudatus'


#.v stands for venous phenotype
group1.v <- as.data.frame(strsplit(group1.v, "\n"))
group1.v$Species <- sub(" ", "_", group1.v$c..Colobus.angolensis....Pygathrix.nemaeus....Nasalis.larvatus...)
group1.v$Species <- sub(" ", "", group1.v$Species)
group1.v <- group1.v %>% select(Species)
#retain only species present in phenotypes data frame
species1.v <- df %>% filter(df$Species %in% group1.v$Species)
species1.v <- species1.v[, c(1, 5:ncol(species1.v))]

#calculate mean for each mirna family (percentage occurrence)
##
means.v <- colMeans(species1.v[, 2:ncol(species1.v)], na.rm = T) # calculate means, results will be in new df)
means_row.v <- data.frame(occurrence.perc = means.v) #
means_row.v <- means_row.v %>% add_row(.before = 1) #adding a row to match the num of cols in species1.v 
rownames(means_row.v) <- colnames(species1.v) #match row and col names of the two dfs before merging 
means_row.v[1, 1]<- "occ_perc"
species1.v <- rbind(species1.v, t(means_row.v)) # add the new means_row.v into the species1.v df

#change species col into rownames
species1.v <- as.data.frame(species1.v)
rownames(species1.v) <- species1.v$Species
species1.v <- species1.v[, c(2:ncol(species1.v))]
##################################################################################MAIN STEP
#filter only mirnas (cols), which have occurrence == 100% across all species
species1.v.100mirna <- as.data.frame(t(species1.v))
species1.v.100mirna <- species1.v.100mirna %>% filter(occ_perc == 1)
species1.v_list100 <- as.list(rownames(species1.v.100mirna)) # save names of mirnas occurring in species 1 (with occurrence100%), as a list
# ↑↑↑↑↑ use space above to get different occurrence percentages i.e, >=90% ↑↑↑↑↑
#################################################################################MAIN STEP

group2.v <- 'Catagonus wagneri
Sus scrofa
Tragulus javanicus
Giraffa camelopardalis
Okapia johnstoni
Rangifer tarandus
Odocoileus virginianus
Elaphurus davidianus
Cervus elaphus
Hemitragus hylocrius
Orvx dammah
Hippotragus niger
Hersact wesuinus
Tragelaphus eurycerus
Bubalus buballs
Bos taurus
Moschus moschiferus
Antilocapra americana
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
Equus asinus'

#.m stands for membrane phenotype
group2.v <- as.data.frame(strsplit(group2.v, "\n"))
group2.v$Species <- sub(" ", "_", group2.v$c..Catagonus.wagneri....Sus.scrofa....Tragulus.javanicus....Giraffa.camelopardalis...)
group2.v$Species <- sub(" ", "", group2.v$Species)
group2.v <- group2.v %>% select(Species)
species2.v <- df %>% filter(df$Species %in% group2.v$Species)
species2.v <- species2.v[, c(1, 5:ncol(species2.v))]

#calculate mean for each mirna family (percentage occurrence)
means.v <- colMeans(species2.v[, 2:ncol(species2.v)], na.rm = T) # calculate means, results will be in new df)
means_row.v <- data.frame(occurrence.perc = means.v) #
means_row.v <- means_row.v %>% add_row(.before = 1) #adding a row to match the num of cols in species2.v 
rownames(means_row.v) <- colnames(species2.v) #match row and col names of the two dfs before merging 
means_row.v[1, 1]<- "occ_perc"
species2.v <- rbind(species2.v, t(means_row.v)) # add the new means_row.v into the species2.v df

#change species col into rownames
species2.v <- as.data.frame(species2.v)
rownames(species2.v) <- species2.v$Species
species2.v <- species2.v[, c(2:ncol(species2.v))]

#filter only mirnas (cols), which have occurrence == 100% across all species
species2.v.100mirna <- as.data.frame(t(species2.v))
species2.v.100mirna <- species2.v.100mirna %>% filter(occ_perc == 1)
species2.v_list100 <- as.list(rownames(species2.v.100mirna))


#[2]
#now obtain: mirnas in 1 but not in 2, and mirnas in 2, but not in 1
## unique.species1.100occ <- mirnas unique to species 1 with 100% occurrence across all species in group 1
## unique.species2.100occ <- mirnas unique to species 1 with 100% occurrence across all species in group 2
## both mirnas present in group 1 and group 2

unique.species1.100occ.v <- setdiff(species1.v_list100, species2.v_list100)
unique.species2.100occ.v  <- setdiff(species2.v_list100, species1.v_list100)
similar.species1.2.100occ.v <- intersect(species2.v_list100, species1.v_list100)


#save the unique lists
#some modifications before saving 
unique.species1.100occ.v <- t(as.data.frame(unique.species1.100occ.v))
unique.species2.100occ.v <- t(as.data.frame(unique.species2.100occ.v))
similar.species1.2.100occ.v <- t(as.data.frame(similar.species1.2.100occ.v))

#save!
'write.csv(unique.species1.100occ.v, file = "C:/Users/jaxk2/Downloads/unique.species1.100occ.v", row.names = F , quote = F)
write.csv(unique.species2.100occ.v, file = "C:/Users/jaxk2/Downloads/unique.species2.100occ.v", row.names = F, quote = F )
write.csv(similar.species1.2.100occ.v, file = "C:/Users/jaxk2/Downloads/similar.species1.2.100occ.v", row.names = F, quote = F)'




#end result: three lists of miRNAs (unique to gorup1, unique to group2, shared between the two gorups)


#[C] there are miRNA lists for each group, so how do the alignments between these groups look like?
# this code cobines each group mirnas with each other, creaitng all possible pairings between the two 
# for each pair, mature miRNA alignement values (dissimilarity from part 1) are obtained
# 3 alignment values are then compared with Krusk-Wallis test, and Dunn test (if Krusk Wallis is significant)
#####################################################################################################################################################
grp1 <- read.csv("C:/Users/jaxk2/Downloads/unique.species1.100occ.v")
grp2 <- read.csv("C:/Users/jaxk2/Downloads/unique.species2.100occ.v")
shared <-read.csv("C:/Users/jaxk2/Downloads/similar.species1.2.100occ.v")
#read in the pairwise mirna co-occurrence distance and dissimilarity average lists 
merged.df.pt1 <- read.csv("C:/Users/jaxk2/Downloads/merged.df.pt1.csv", header = T)
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
g1.vs.g2_merged.df.pt1 <- na.omit(g1.vs.g2_merged.df.pt1) #some mirnas are not present in the abs/pres matrix, hence to get rid of those do na.omit
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#[2] group1 mutually exclusive vs shared mirnas
####################################################################################################################################################################
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#[3] group2 mutually exclusive vs shared mirnas
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#[4] plot the dissimilarity values for the 3 different groups
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#see dissimilarity distribution for each group with 1s
##first group_by the categorical variables (i.e., the 3 groups) 
dissimilarity_3groups <- dissimilarity_3groups %>%
  group_by(label)


ll <- dissimilarity_3groups %>%
  ggplot(aes(x = dissimilarity_average)) +
  geom_histogram(aes(fill = label), bins = 25) +
  facet_wrap(~ label) +
  labs(title = "miRNA Dissimilarity distribution (Venous Pattern)", 
       x = "miRNA sequence Dissimilarity", 
       fill = 'Groups') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) #plot title in the center


#ggsave("C:/Users/jaxk2/Downloads/Venous_dissimilarity_distribution_supp.png", plot = ll, width = 6, height = 4, units = "in", dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~ for proportion of zeros in each group  
#proportion of 0s in each group
dissimilarity_3groups <- dissimilarity_3groups <- dissimilarity_3groups %>%
  mutate(is_zero = dissimilarity_average == 0)

proportions <-  dissimilarity_3groups %>%
  group_by(label, is_zero) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

kk <- ggplot(proportions, aes(x = label, y = proportion, fill = as.factor(is_zero))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("blue", "orange"), labels = c("Identical", "Non-Identical")) +
  labs(x = "Comparison Groups", y = "Proportion", fill = " ", 
       title = 'Proportion of Identical and non-Identical Sequences') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

#ggsave("C:/Users/jaxk2/Downloads/Venous_proportions_MAIN.png", plot = kk, width = 6, height = 4, units = "in", dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#[5] kruskal and DUNN's test 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
krusk <- kruskal.test(dissimilarity_average ~ label, data = dissimilarity_3groups)
dunn <- dunnTest(dissimilarity_average ~ label, data = dissimilarity_3groups, method = "bonferroni")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#[6] proportions test
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
g1g2_vs_g1shared_test <- prop.test(c(29, 264), c(13+29, 142+264)) #X-squared = 0.12348, df = 1, p-value = 0.7253
g1g2_vs_g2shared_test <- prop.test(c(29, 89), c(13+29, 62+89)) #X-squared = 1.0195, df = 1, p-value = 0.3126
g1shared_vs_g2shared_test <- prop.test(c(264, 89), c(142+264, 62+89)) #X-squared = 1.503, df = 1, p-value = 0.2202

#none is significant 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

