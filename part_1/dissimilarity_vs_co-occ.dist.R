library(spaa) # packageVersion('spaa') '0.2.2'
library(tidyverse) # packageVersion('tidyverse') '2.0.0'
library(dplyr) # packageVersion('dplyr') ‘1.1.3’
library(vegan) # packageVersion('vegan') ‘2.6.4’
library(ggplot2) # packageVersion('ggplot2') ‘3.4.4’
library(readxl) #packageVersion('readxl') ‘1.4.3’
library(ISLR) #packageVersion('ISLR') ‘1.4’
library(ggpubr) # packageVersion('ggpubr') ‘0.6.0’
library(FSA) # for DUNN test # packageVersion('FSA') ‘0.9.5’
library(ggsignif) # packageVersion('ggsignif') ‘0.6.4’
#every section separated by 5 lines and numbered in []





#[1] MIRNA STRING LABELS MODIFICATION
setwd("C:/Users/jaxk2")
getwd()
# code run with all different alignments, input file commented in line 21
alignment_wd4 <- read.table('C:/Users/jaxk2/Downloads/all_fasta_1seq_alignment') #alignment4_mirnas_of_interest_wd4 #all_fasta_1seq_alignment #alignment2 
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
#intersect (no need as the alignment was conducted using the updated mirnas of interest to seqkit the fasta file. so that the fasta only had mirnas from the co-occurrence matrix)
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

#now calculate the miRNA co-occurrence distance for the updated_df.csv
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

#average the alignments, as some sequences aligned multiple times
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
## values of 1 for dissimilarity are manually inserted as 30 random 30 non-aligned sequences were non-identical when manually aligned
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




#create a column with dissimilarity group ranges: 
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
sum(is.na(merged.df.na$dissimilarity_average)) #154755 is.na (i.e., not aligned)


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
  select(merged_mirna, mirna1.x, mirna2.x,distance_score, dissimilarity_average, avg.diss.groups, avg.distance.groups)

merged.df.na.omit <- na.omit(merged.df.na)
anyNA(merged.df.na.omit$avg.distance.groups) #FALSE for all 3 groups------------#check





#[7]
#Viewing some plots: distributions, boxplots, integrating p values in plots

#-----------------------------DISTRIBUTIONS-------------------------------------start

## DISSIMILARITY DISTRIBUTION
#with 1s added for non-aligned
s1 <- ggplot(merged.df, aes(x = dissimilarity_average)) +
  geom_histogram(fill = 'orange' ) +
  labs(x = "miRNA Sequence Dissimilarity", 
       y = "Frequency", 
       title = "Distribution miRNA Sequence Dissimilarity") +
  theme_minimal(12) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect()) +
  theme(plot.title = element_text(hjust = 0.5))

# save high quality figures
ggsave("C:/Users/jaxk2/Downloads/dissimilarity_distribution_supp.png", plot = s1, width = 6, height = 4, units = "in", dpi = 300)


#distribution of seq dissimilarity score for aligned pairs (pairs in both alignments and miRNA co-occurrence df)
f2_dd <- ggplot(merged.df.na.omit, aes(x = dissimilarity_average)) +
  geom_histogram(fill = 'blue') +
  labs(x = "miRNA Sequence Dissimilarity", 
       y = "Frequency", 
       title = "Distribution of miRNA Sequence Dissimilarity") +
  theme_minimal()  +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect()) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("C:/Users/jaxk2/Downloads/dissimilarity_distribution_MAIN.png", plot = f2_dd, width = 6, height = 4, units = "in", dpi = 300)



## DISTANCE SCORE DISTRIBUTIONS
#with all miRNA pairs  
s2 <- ggplot(merged.df, aes(x = distance_score)) +
  geom_histogram(fill = 'orange') +
  labs(x = "miRNA Co-Occurrence Distance", 
       y = "Frequency", 
       title = "Distribution of miRNA Co-Occurrence Distance") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect())  +
  theme(plot.title = element_text(hjust = 0.5)) #plot title in the center

ggsave("C:/Users/jaxk2/Downloads/distance_distribution_supp.png", plot = s2, width = 6, height = 4, units = "in", dpi = 300)


#with aligned miRNA pairs only
f2_cd <- ggplot(merged.df.na.omit, aes(x = distance_score)) +
  geom_histogram(fill = 'blue') +
  labs(x = "miRNA Co-Occurrence Distance", 
       y = "Frequency", 
       title = "Distribution of miRNA Co-Occurrence Distance") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect()) +
  theme(plot.title = element_text(hjust = 0.5)) #plot title in the center

ggsave("C:/Users/jaxk2/Downloads/distance_distribution_MAIN.png", plot = f2_cd, width = 6, height = 4, units = "in", dpi = 300)


#for all distance scores for all pairs, before merging the alignment pairs
s2.5 <- ggplot(mirna_dist, aes(x = distance_score)) +
  geom_histogram() +
  labs(x = "miRNA Co-Occurrence Distance", 
       y = "Frequency", 
       title = "Distribution of miRNA Co-Occurrence Distance") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  theme(panel.background = element_rect()) +
  theme(plot.title = element_text(hjust = 0.5)) #plot title in the center

ggsave("C:/Users/jaxk2/Downloads/distance_distribution.all_supp.png", plot = s2.5, width = 6, height = 4, units = "in", dpi = 300)

#-----------------------------DISTRIBUTIONS-------------------------------------end


#-------------------------BOX-PLOTS/REGRESSIONS---------------------------------start
#to plot relationship between seq dissimilarity and miRNA absence/presence. 

#REGRESSION LINE 
#only with aligned 
f <- ggplot(merged.df.na.omit, aes(x = distance_score, y = dissimilarity_average)) +
  geom_point(alpha=10) +
  labs(x = "miRNA Co-occurence Distance", 
       y = "miRNA Sequence Dissimilarity", 
       title = "miRNA Co-Occurrence Distance vs. miRNA Sequence Dissimilarity") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) #plot title in the center
# save at higher quality
ggsave("C:/Users/jaxk2/Downloads/correlation_MAIN.png", plot = f, width = 6, height = 4, units = "in", dpi = 300)

#with non-aligned
s3 <- ggplot(merged.df, aes(x = distance_score, y = dissimilarity_average)) +
  geom_point(alpha=0.9) +
  labs(x = "miRNA Co-occurence Distance", 
       y = "miRNA Sequence Dissimilarity", 
       title = "miRNA Co-ccurrence Distance vs. miRNA Sequence Dissimilarity") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) #plot title in the center

ggsave("C:/Users/jaxk2/Downloads/correlation_supp.png", plot = s3, width = 6, height = 4, units = "in", dpi = 300)


#BOXPLOTS

#x = dissimilarity groups, y = distance scores
# used in main thesis

#aligned + non aligned
s_boxplot <- ggplot(merged.df) +
  aes(x = avg.diss.groups, y = distance_score, fill = avg.diss.groups) +
  geom_boxplot() +
  labs(x = 'miRNA Sequence Dissimilarity', #average mature miRNA family pairwise dissimilarity score groups 
       y = 'miRNA Co-Occurrence Distance', 
       title = 'miRNA Sequence Dissimilarity groups vs
       miRNA Co-Occurrence Distance', 
       fill = 'miRNA sequence
Dissimilarity groups') +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5)) + #plot title in the center
  scale_fill_brewer(palette = "Set1") #color-blind friendly palette


#only aligned
p <- ggplot(merged.df.na.omit) +
  aes(x = avg.diss.groups, y = distance_score, fill = avg.diss.groups) +
  geom_boxplot() +
  labs(x = 'miRNA Sequence Dissimilarity', #average mature miRNA family pairwise dissimilarity score groups 
       y = 'miRNA Co-Occurrence Distance', 
       title = 'miRNA Sequence Dissimilarity groups vs
       miRNA Co-Occurrence Distance', 
       fill = 'miRNA sequence
Dissimilarity groups') +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5))  +
  scale_fill_brewer(palette = "Set1") + #color-blind friendly palette
  theme_bw()
#-------------------------BOX-PLOTS/REGRESSIONS---------------------------------end





#[8] krusk and Dunn's tests
#useful for results section
doBy::summary_by(distance_score ~ avg.diss.groups,
                 data = merged.df, 
                 FUN = median)

# Kruskal Wallis adn Dunn's test with non aligned sequences included				 
kwt <- kruskal.test(distance_score ~ avg.diss.groups, data = merged.df)
dunn.t <-  dunnTest(distance_score ~ avg.diss.groups, data = merged.df, method = "bonferroni")

# Kruskal Wallis adn Dunn's test with only aligned sequences(n main thesis)
kwtkwt.na.omit <- kruskal.test(distance_score ~ avg.diss.groups, data = merged.df.na.omit)
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



#add significant p-values from Dunn's test to boxplot
#use previous ggplot to add significance values; from line 360 

pp <- p + geom_signif(
  comparisons = list(
    c("0%", "15-20%"),
    c("0%", "5-10%")
  ),
  map_signif_level = TRUE,
  annotations = c("**", "*"),
  y_position = c(max(merged.df.na.omit$distance_score) + 0.1, max(merged.df.na.omit$distance_score) + 0.2),
  tip_length = 0.03
)

ggsave("C:/Users/jaxk2/Downloads/Boxplot_sig_MAIN.png", plot = pp, width = 6, height = 4, units = "in", dpi = 300)




#[10] is there a difference in the number of sequences with 0 dissimilarity score between co-occurring and not co-occurring sequences? --> carry out a prop.test
# co-occurring sequences (distance scores == 0)  ;  not co-occurring sequences (distance score == 1)
# conduct analysis for only aligned sequences and for all sequences (aligned+non-aligned) 

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

data:  c(901, 5582) out of c(185 + 901, 5582 + 2105)
X-squared = 52.299, df = 1, p-value = 4.767e-13
alternative hypothesis: two.sided
95 percent confidence interval:
 0.07848307 0.12849501
sample estimates:
   prop 1    prop 2 
0.8296501 0.7261611 '

#plot the prop data
#create a column with proportion values for each group
prop_cont$proportions <- c(185/(185+901), 2105/(2105+5582), 901/(185+901), 5582/(2105+5582))

prop_hist <- ggplot(prop_cont, aes(x = Var1, y = proportions, fill = as.factor(Var2))) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c("blue", "orange"), labels = c("Identical \nsequences", "\nNon-Identical \nsequences")) +
  labs(x = "Groups", y = "Proportion", fill = " ",
       title = "Proportion of Similar Sequences in 
Co-Occurring vs non-Co-Occurring miRNA Pairs") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons = list(c("co-occurring", "not co-occurring")),
              map_signif_level = TRUE,
              annotations = "***",
              y_position = max(prop_cont$proportions) + 0.2,
              tip_length = 0.03)

ggsave("C:/Users/jaxk2/Downloads/Proportion_hist_MAIN.png", plot = prop_hist, width = 5.8, height = 5, units = "in", dpi = 300)






# save a list of non-aligned miRNA pairs
## randomly align 30 of these pairs on https://www.ebi.ac.uk/jdispatcher/psa with local alignment 

difference <- setdiff(mirna_dist$merged_mirna, alignment_wd4$merged_mirna)
difference <- as.data.frame(difference)

difference <- setNames(difference, 'unaligned_pairs')
difference <- difference %>% 
  mutate(mirna1 = sub("-.*", "", unaligned_pairs), 
         mirna2 = sub(".*-", "", unaligned_pairs))

#generate 30 random numbers
random_nums <- as.array(sample(1:154755, 30))
#extract the numbers of row matching with the random numbers form the non-aligned mirnas df (difference)
random_difference <- difference[random_nums, ]
#save
'write.csv(random_difference,file = "C:/Users/jaxk2/Downloads/30_random_non-aligned_1seq.csv", row.names = FALSE, col.names = FALSE)'


#Save df with miRNA co-occurrence distance scores and alignment dissimialrity scores
# To be used in analysis for part 2 ...
'write.csv(merged.df, file = "C:/Users/jaxk2/Downloads/merged.df.pt1_1seq", col.names = T)'
