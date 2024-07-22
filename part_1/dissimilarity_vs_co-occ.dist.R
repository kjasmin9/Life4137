library(spaa)
library(tidyverse)
library(dplyr)
library(vegan)
library(ggplot2)
library(readxl)
library(ISLR)
library(ggpubr)
library(FSA) # for DUNN test
library(ggsignif)
#every section separated by 5 lines and numbered in []





#[1] MIRNA STRING LABELS MODIFICATION
#alignment3_wd4 <- read.table('C:/Users/jaxk2/Downloads/alignment3_wd4') #older alignement, has way less alignemnts
alignment_wd4 <- read.table('C:/Users/jaxk2/Downloads/all_fasta_1seq_alignment') #alignment4_mirnas_of_interest_wd4 #all_fasta_1seq_alignment alignment_100_not100seqs
colnames(alignment_wd4) <- c('column1', 'column2', 'column3')

#match miRNA names strings to the updated_df
alignment_wd4$mirna1 <- sub(".*?-", "", alignment_wd4$column1) #delete everything before the first '-' (inclusive)
alignment_wd4$mirna1 <- sub("^(([^-]*-){1}[^-]*)-.*", "\\1", alignment_wd4$mirna1) #delete everything after the 2nd '-' (inclusive)
alignment_wd4$mirna1 <- gsub("-", "", alignment_wd4$mirna1) #remove the '-' 
alignment_wd4$mirna1 <- sub("_.*", "", alignment_wd4$mirna1) #delete everything after the '_' (inclusive) 

alignment_wd4$mirna2 <- sub(".*?-", "", alignment_wd4$column2) #delete everything before the first '-' (inclusive)
alignment_wd4$mirna2 <- sub("^(([^-]*-){1}[^-]*)-.*", "\\1", alignment_wd4$mirna2) #delete everything after the 2nd '-' (inclusive)
alignment_wd4$mirna2 <- gsub("-", "", alignment_wd4$mirna2) #remove the '-' 
alignment_wd4$mirna2 <- sub("_.*", "", alignment_wd4$mirna2) #delete everything after the '_' (inclusive) 

#n_distinct(alignment_wd4$mirna1) # there are 431 mirna families in the alignment #1_seq 430 mirnas
#unique(alignment_wd4$mirna1) #prints all the unique mirnas (430)





#[2] ALIGNMENT ONLY TO CONTAIN mirnas IN UPDATED_DF
#intersect (no need because the aligenment was consucted using the updated mirnas of interest to seqkit the fasta file. so tha tthe fasta only had mirnas form the prs/abs matrix)
mirnas_of_interest <- read.csv('C:/Users/jaxk2/OneDrive/Documents/Downloads/updated_df1.csv', header = F )
mirnas_of_interest <- mirnas_of_interest[1, c(5:ncol(mirnas_of_interest))]
mirnas_of_interest <- t(mirnas_of_interest)

#------------------------------------------------------------------------------#
all(mirnas_of_interest[,1]%in%alignment_wd4$mirna1) #FALSE                      #check
all(mirnas_of_interest[,1]%in%alignment_wd4$mirna2) #FALSE                      #check  
all(alignment_wd4$mirna1%in%mirnas_of_interest[,1]) #TRUE                       #check
all(alignment_wd4$mirna2%in%mirnas_of_interest[,1]) #TRUE                       #check  
all(alignment_wd4$mirna1%in%alignment_wd4$mirna2) #TRUE                         #check
#conclusions: 
#          1)not all mirnas form the updated_df.csv are present in the alignment
#          2)code worked: alignment has only mirnas from updated_df.csv 
#          3)in the alignment, both mirna1 and 2 have the same mirnas set
#------------------------------------------------------------------------------#

setdiff(mirnas_of_interest, alignment_wd4$mirna1) # to see the miRNA families that are not in the alignement file, do length at beginning to get the numeber (52)

#matching updated_df.csv mirnas to the set of mirnas present in the alignemnt 
miRNA_data <- read.csv('C:/Users/jaxk2/OneDrive/Documents/Downloads/updated_df1.csv', header = F) #read in the absence/presence df
miRNA_data <- miRNA_data[, c(5:ncol(miRNA_data))]
miRNA_data<- t(miRNA_data)
miRNA_data <- as.data.frame(miRNA_data)

#filter, so that the updated.csv has only mirnas present in the alignment 
miRNA_data1 <- miRNA_data %>% filter(miRNA_data[, 1]%in%alignment_wd4$mirna1) 

#------------------------------------------------------------------------------##check
all(miRNA_data1[,1]%in%alignment_wd4$mirna1) #TRUE                              #check
all(alignment_wd4$mirna1%in%miRNA_data1[,1]) #TRUE                              #check
#conclusions 
#            alignemnt & updated.csv have same set of miRNAs
#------------------------------------------------------------------------------##check





#[3] DISTANCE calculation for miRNA co-occurrence
#set row names for the miRNA-data1
rownames(miRNA_data1) <- miRNA_data1[, 1]
miRNA_data1 <- miRNA_data1[, c(2:ncol(miRNA_data1))]

#now calculate the miRNA co-occurence distance for the updated_df.csv
mirna_dist <- dist(miRNA_data1, method = 'binary', na.omit(TRUE)) 
dim(as.matrix(mirna_dist)) #430 x 430 
anyNA(mirna_dist) #check, should be FALSE 





#[4] MERGING THE sequence alignment(alignment_wd4) with the mirna co-occurrence distance matrix(mirna_dist)
#convert dist matrix into list format
mirna_dist <- dist2list(mirna_dist)
#colnames
mirna_dist <- setNames(mirna_dist, c('mirna1', 'mirna2', 'distance_score'))
alignment_wd4 <- setNames(alignment_wd4, c('ogmirna1','ogmirna2', 'seq_similarity','mirna1', 'mirna2' ))
alignment_wd4$dissimilarity_perc <- (100 - alignment_wd4$seq_similarity)/100

#create a column of reference for both miRNA_dist and alignemnt_wd4
merged_mirna <- mirna_dist %>% unite(merged_mirna, mirna1, mirna2, sep = '-')
mirna_dist$merged_mirna <- merged_mirna$merged_mirna #addcol
merged_mirna2 <- alignment_wd4 %>% unite(merged_mirna, mirna1, mirna2, sep = '-')
alignment_wd4$merged_mirna <- merged_mirna2$merged_mirna #addcol

n_distinct(alignment_wd4$merged_mirna) #check #[1] 30252

#average the alignemnts
avg.alignment_wd4 <-  alignment_wd4 %>%
  group_by(merged_mirna) %>%
  summarise(dissimilarity_average = mean(dissimilarity_perc, na.rm = T)
  )
#separate the mirna pairs into two cols. 
avg.alignment_wd4 <- avg.alignment_wd4 %>%
  mutate(mirna1 = sub("-.*", "", merged_mirna), 
         mirna2 = sub(".*-", "", merged_mirna))


#merge the two: co-occurrence distance and alignment similarity scores
## keep 2 merged.dfs one with NA values and one where the NA values are replaced by 1s (in avg.dissimilarity)
## values of 1 for dissimilarity are amually inserted as 30 random 30 non-aligned sequences were non-identical when manually aligned
### use mirna_dist and avg.alignemnt_wd4

merged.df.na <- merge(mirna_dist, avg.alignment_wd4, by = "merged_mirna", all.x = T)
sum(!is.na(merged.df.na$dissimilarity_average)) # out of 184,900 mirna family pairs:  154755 not aligned (NA), 30,145 aligned
merged.df <- merge(mirna_dist, avg.alignment_wd4, by = "merged_mirna", all.x = T)
sum(is.na(merged.df$dissimilarity_average)) # 171687                            #check
merged.df$dissimilarity_average[is.na(merged.df$dissimilarity_average)] <- 1
#save to use in ada
'write.csv(merged.df, file = "C:/Users/jaxk2/Downloads/merged.df_1seq", col.names = T, row.names = F)
write.csv(merged.df.na, file = "C:/Users/jaxk2/Downloads/merged.df.na_1seq", col.names = T, row.names = F)'





#[5] on ada: list2dist.R, sbatch_4_r.sh, mantel_test.R 





#[6] Grouping sequence dissimilarity and co-occurrence distance scores into ranges
##use merged.df and merged.df.na

#to group you first need to see the ranges of the dissimilarity values 
  
range(merged.df$dissimilarity_average) #0.0 1
range(na.omit(merged.df.na$dissimilarity_average)) # 0.0 0.2
summary(merged.df$dissimilarity_average)
summary(na.omit(merged.df.na$dissimilarity_average))
summary(merged.df$distance_score) 
summary(merged.df.na$distance_score)



  
#create column with dissimilarity group ranges: 
merged.df <- merged.df %>% 
  mutate(avg.diss.groups = case_when(dissimilarity_average == 0  ~ "0%", 
                                     dissimilarity_average > 0 & dissimilarity_average <= 0.05 ~ "0-5%", 
                                     dissimilarity_average > 0.05 & dissimilarity_average <= 0.1 ~ "5-10%",
                                     dissimilarity_average > 0.1 & dissimilarity_average <= 0.15 ~ "10-15%",
                                     dissimilarity_average > 0.15 & dissimilarity_average <= 0.2 ~ "15-20%",
                                     dissimilarity_average > 0.99 ~ "100%",))



merged.df.na <- merged.df.na %>% 
  mutate(avg.diss.groups = case_when(dissimilarity_average == 0  ~ "0%", 
                                     dissimilarity_average > 0 & dissimilarity_average <= 0.05 ~ "0-5%", 
                                     dissimilarity_average > 0.05 & dissimilarity_average <= 0.1 ~ "5-10%",
                                     dissimilarity_average > 0.1 & dissimilarity_average <= 0.15 ~ "10-15%",
                                     dissimilarity_average > 0.15 & dissimilarity_average <= 0.2 ~ "15-20%",))


sum(!is.na(merged.df.na$dissimilarity_average)) #30,145 !is.na (i.e. aligned) 154755 is.na (i.e., not aligned)
sum(!is.na(merged.df.na$dissimilarity_average)) #154755 is.na (i.e., not aligned)


#repeat with distance scores
merged.df <- merged.df %>% 
  mutate(avg.distance.groups = case_when(distance_score == 0  ~ "0%", 
                                         distance_score > 0 & distance_score <= 0.2 ~ "0-20%", 
                                         distance_score > 0.2 & distance_score <= 0.4 ~ "20-40%",
                                         distance_score > 0.4 & distance_score <= 0.6 ~ "40-60%",
                                         distance_score > 0.6 & distance_score <= 0.8 ~ "60-80%",
                                         distance_score > 0.8 & distance_score < 1 ~ "80-99%", 
                                         distance_score >= 1  ~ "100%")) 


merged.df.na <- merged.df.na %>% 
  mutate(avg.distance.groups = case_when(distance_score == 0  ~ "0%", 
                                         distance_score > 0 & distance_score <= 0.2 ~ "0-20%", 
                                         distance_score > 0.2 & distance_score <= 0.4 ~ "20-40%",
                                         distance_score > 0.4 & distance_score <= 0.6 ~ "40-60%",
                                         distance_score > 0.6 & distance_score <= 0.8 ~ "60-80%",
                                         distance_score > 0.8 & distance_score < 1 ~ "80-99%", 
                                         distance_score >= 1  ~ "100%")) 



#organise the levels in the new cols (needed for uniformity in boxplots)

merged.df <- merged.df %>%
  mutate(avg.diss.groups = factor(avg.diss.groups, levels = c("0%", "0-5%", "5-10%", "10-15%", "15-20%", "100%")))

merged.df <- merged.df %>%
  mutate(avg.distance.groups = factor(avg.distance.groups, levels = c("0%", "0-20%", "20-40%", "40-60%", "60-80%", "80-99%", "100%")))

merged.df.na <- merged.df.na %>%
  mutate(avg.diss.groups = factor(avg.diss.groups, levels = c("0%", "0-5%", "5-10%", "10-15%", "15-20%")))

merged.df.na <- merged.df.na %>%
  mutate(avg.distance.groups = factor(avg.distance.groups, levels = c("0%", "0-20%", "20-40%", "40-60%", "60-80%", "80-99%", "100%")))

merged.df.na.omit <- merged.df.na %>%
  select(merged_mirna, mirna1.x, mirna2.x,distance_score, dissimilarity_average, avg.diss.groups, avg.diss.groups2, avg.distance.groups, avg.distance.groups2 )

merged.df.na.omit <- na.omit(merged.df.na)
anyNA(merged.df.na.omit$avg.distance.groups) #FALSE for all 3 groups------------#check





#[7]
#Viewing some plots: distributions, boxplots, integrating p values in plots

#-----------------------------DISTRIBUTIONS-------------------------------------start

## DISSIMILARITY DISTRIBUTION
#with 1s added for non-aligned
ggplot(merged.df, aes(x = dissimilarity_average)) +
  geom_histogram(fill = 'blue' ) +
  labs(x = "Pairwise mature miRNA sequence dissimilarity", 
       y = "Frequency", 
       title = "Distribution of Pairwise Mature miRNA Sequence Dissimilarity") +
  theme_minimal(12) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect())
  
#distribution of seq dissimilarity score for aligned pairs (pairs in both alignments and miRNA co-occurrence df)
ggplot(merged.df.na.omit, aes(x = dissimilarity_average)) +
  geom_histogram(fill = 'orange') +
  labs(x = "Pairwise mature miRNA sequence dissimilarity", 
       y = "Frequency", 
       title = "Distribution of Pairwise Mature miRNA Sequence Dissimilarity") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect())
  

## DISTANCE SCORE DISTRIBUTIONS
#with all miRNA pairs  
ggplot(merged.df, aes(x = distance_score)) +
  geom_histogram(fill = 'blue') +
  labs(x = "Distribution of miRNA co-occurrence Distance Scores", 
       y = "Frequency", 
       title = "miRNA co-occurrence Distance Scores") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect())

#with aligned miRNA pairs only
ggplot(merged.df.na.omit, aes(x = distance_score)) +
  geom_histogram(fill = 'orange') +
  labs(x = "Distribution of miRNA co-occurrence Distance Scores", 
       y = "Frequency", 
       title = "miRNA co-occurrence Distance Scores") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect())

#for all distance scores for all pairs, before merging the alignment pairs
ggplot(mirna_dist, aes(x = distance_score)) +
  geom_histogram() +
  labs(x = "Distribution of miRNA co-occurrence Distance Scores", 
       y = "Frequency", 
       title = "miRNA co-occurrence Distance Scores") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect())

#-----------------------------DISTRIBUTIONS-------------------------------------end


#-------------------------BOX-PLOTS/REGRESSIONS---------------------------------start
#to plot relationship between seq dissimilarity and miRNA absence/presence. 

#REGRESSION LINE 
#only with aligend 
ggplot(merged.df.na.omit, aes(x = distance_score, y = dissimilarity_average)) +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm",formula = y~poly(x,4)) +
  labs(x = "miRNA Co-occurence Distance Score", 
       y = "Pairwise mature miRNA sequence dissimilarity", 
       title = "Relationship between miRNA Co-occurence Distance Score & 
       Pairwise mature miRNA sequence dissimilarity") +
  theme_classic()

#with non-aligned
ggplot(merged.df, aes(x = distance_score, y = dissimilarity_average)) +
  geom_point(alpha=0.9) +
  stat_smooth(method="lm",formula = y~poly(x,4))+
  labs(x = "miRNA Co-occurence Distance Score", 
       y = "Pairwise mature miRNA sequence dissimilarity", 
       title = "Relationship between miRNA Co-occurence Distance Score & 
       Pairwise mature miRNA sequence dissimilarity") +
  theme_classic()


#BOXPLOTS


#aligned + non aligned
##avg.***groups
#x = dissimilarity groups, y = distance scores
ggplot(merged.df) +#input the 3 different: merged.df, merged.df.na, merged.df.na.omit
  aes(x = avg.diss.groups, y = distance_score, fill = avg.diss.groups) +
  geom_boxplot() +
  labs(x = 'Mature Sequence dissimilarity between pairs of miRNAs ', #average mature miRNA family pairwise dissimilarity score groups 
       y = 'Level of co-occurrence between pair of miRNA families', 
       title = 'Pairwise mature miRNA sequence dissimilarity groups &
       Level of co-occurrence between pair of miRNA families boxplots') +
  theme(legend.position = "none")

#x = distance_groups y = sissimilarity scores
ggplot(merged.df) +#input teh 3 different: merged.df, merged.df.na, merged.df.na.omit #na.omit yield interesting results
  aes(x = avg.distance.groups, y = dissimilarity_average, fill = avg.distance.groups) +
  geom_boxplot() +
  labs(x = 'co-occurrence distance ', #average mature miRNA family pairwise dissimilarity score groups 
       y = 'dissimilarity between mature sequences', 
       title = 'Pairwise mature miRNA sequence dissimilarity groups &
       Level of co-occurrence between pair of miRNA families boxplots') +
  theme(legend.position = "none")



#only aligned
#x = dissimilarity groups, y = distance scores
ggplot(merged.df.na.omit) +#input the 3 different: merged.df, merged.df.na, merged.df.na.omit
  aes(x = avg.diss.groups, y = distance_score, fill = avg.diss.groups) +
  geom_boxplot() +
  labs(x = 'Mature Sequence dissimilarity between pairs of miRNAs ', #average mature miRNA family pairwise dissimilarity score groups 
       y = 'Level of co-occurrence between pair of miRNA families', 
       title = 'Pairwise mature miRNA sequence dissimilarity groups &
       Level of co-occurrence between pair of miRNA families boxplots') +
  theme(legend.position = "left")


#x = distance_groups y = dissimilarity scores
ggplot(merged.df.na.omit) +#input teh 3 different: merged.df, merged.df.na, merged.df.na.omit #na.omit yield interesting results
  aes(x = avg.distance.groups, y = dissimilarity_average, fill = avg.distance.groups) +
  geom_boxplot() +
  labs(x = 'co-occurrence distance ', #average mature miRNA family pairwise dissimilarity score groups 
       y = 'dissimilarity between mature sequences', 
       title = 'Pairwise mature miRNA sequence dissimilarity groups &
       Level of co-occurrence between pair of miRNA families boxplots') +
  theme(legend.position = "none")





#[8] krusk and Dunn's tests
#useful for results section
doBy::summary_by(distance_score ~ avg.diss.groups,
                 data = merged.df, 
                 FUN = median)

# Kruskal Wallis adn Dunn's test with non aligned sequences included				 
kwt <- kruskal.test(distance_score ~ avg.diss.groups, data = merged.df)
dunn.t <-  dunnTest(distance_score ~ avg.diss.groups, data = merged.df, method = "bonferroni")

# Kruskal Wallis adn Dunn's test with only aligned sequences 
kwt.na.omit <- kruskal.test(distance_score ~ avg.diss.groups, data = merged.df.na.omit)
dunn.t.na.omit <- dunnTest(distance_score ~ avg.diss.groups, data = merged.df.na.omit, method = "bonferroni")


# arrange Dunn's test's results in a data frame (table to be inserted in figures) 
dunn.t_df <- as.data.frame(dunn.t.na.omit$res)
comparisons <- dunn.t_df$Comparison
p_value <- dunn.t_df$P.adj
dunn.t_df <-  data.frame(
  group1 = sapply(comparisons, function(x) strsplit(x, " - ")[[1]][1]),
  group2 = sapply(comparisons, function(x) strsplit(x, " - ")[[1]][2]),
  p_value = p_value
)

#useful to observe the significance labels to be added to the boxplot
dunn.t_df$label <- ifelse(dunn.t_df$p_value < 0.001, "***",
                              ifelse(dunn.t_df$p_value < 0.01, "**",
                                     ifelse(dunn.t_df$p_value < 0.05, "*", "ns")))
									 


#add significant p-calues from Dunn's test to boxplot
p <- ggplot(merged.df.na.omit) + #store plot as object
  aes(x = avg.diss.groups, y = distance_score, fill = avg.diss.groups) +
  geom_boxplot() +
  labs(
    x = 'Mature Sequence dissimilarity between pairs of miRNAs',
    y = 'Level of co-occurrence between pair of miRNA families',
    title = 'Pairwise mature miRNA Sequence Dissimilarity groups &\nLevel of co-occurrence between pair of miRNA families boxplots',
    fill = 'Dissimilarity groups \nlegend'
  ) +
  theme(legend.position = "right")

p + geom_signif(
  comparisons = list(
    c("0%", "15-20%"),
    c("0%", "5-10%")
  ),
  map_signif_level = TRUE,
  annotations = c("**", "*"),
  y_position = c(max(merged.df.na.omit$distance_score) + 0.1, max(merged.df.na.omit$distance_score) + 0.2),
  tip_length = 0.03
)





#[10] is there a difference in the number of sequences with 0 dissimilarity score between co-occurring and not co-occurring sequences? --> carry out a prop.test
# co-occurring sequences (distance scores == 0)  ;  not co-occurring sequences (distance score == 1)
# conduct analysis for only aligned sequences and for all sequences (aligned+non-aligend) 

prop <- merged.df.na.omit %>%
  select(mirna1.x, mirna2.x, merged_mirna, distance_score, dissimilarity_average)
#retain only values whereby distance is either = to 1 or 0
prop <- prop[prop$distance_score %in% c(0, 1), ]
# add a column to label dissimilarity values by 0 -> 0, if not 0 -> 1 (0 = similar, 1 = not similar
prop <- prop %>%
  mutate(similarity = ifelse(dissimilarity_average == 0, 'similar', 'not_similar'),
         co_occurrence = ifelse(distance_score == 0, 'co-occurring', 'not co-occurring'))

# prepare a contingency table for prop.test 
prop_cont <- table(prop$co_occurrence, prop$similarity)
prop_cont <- as.data.frame(prop_cont[, c(2, 1)])
prop.test(c(901, 5582), c(185+901, 5582+2105))
'2-sample test for equality of proportions with continuity correction

data:  c(901, 5582) out of c(1086, 5582)
X-squared = 971.72, df = 1, p-value < 2.2e-16
alternative hypothesis: two.sided
95 percent confidence interval:
 -0.1932588 -0.1474410
sample estimates:
   prop 1    prop 2 
0.8296501 1.0000000'

#plot the prop data
#create a column with proportion values for each group
prop_cont$proportions <- c(185/(185+901), 2105/(2105+5582), 901/(185+901), 5582/(2105+5582))
ggplot(prop_cont, aes(x = Var1, y = proportions, fill = as.factor(Var2))) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c("blue", "orange"), labels = c("Identical sequences", "Non-Identical sequences")) +
  labs(x = "Groups", y = "Proportion", fill = "Legend") + 
  geom_signif(comparisons = list(c("co-occurring", "not co-occurring")),
              map_signif_level = TRUE,
              annotations = "***",
              y_position = max(prop_cont$proportions) + 0.2,
              tip_length = 0.03)





# save a list of non-aligned miRNA pairs
## randomly align 30 of these pairs on https://www.ebi.ac.uk/jdispatcher/psa with local alignment 

difference <- setdiff(mirna_dist$merged_mirna, alignment_wd4$merged_mirna)
difference <- as.data.frame(difference)

difference <- setNames(difference, 'unaligned_pairs')
difference <- difference %>% 
  mutate(mirna1 = sub("-.*", "", unaligned_pairs), 
         mirna2 = sub(".*-", "", unaligned_pairs))
#save
write.csv(difference,file = "C:/Users/jaxk2/Downloads/non_aligned_pairs_1seq.csv", row.names = FALSE, col.names = FALSE)


#Save df with miRNA co-occurrence distance scores and alignment dissimialrity scores
# To be used in analysis for part 2 ...
write.csv(merged.df, file = "C:/Users/jaxk2/Downloads/merged.df.pt1_1seq", col.names = T)
