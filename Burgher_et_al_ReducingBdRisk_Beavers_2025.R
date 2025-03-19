library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(openxlsx)
library(lme4)
library(plotrix)
detach("package:Rmisc", unload = TRUE)
detach("package:plyr", unload = TRUE)

Burgher_et_al_ReducingBdRiskBeaver_Data <- read_excel("Burgher et al_ReducingBdRiskBeaver_Data.xlsx")

Final_Data <- Burgher_et_al_ReducingBdRiskBeaver_Data

head(Final_Data)

#Want to add 1 in growth column for all rafts and all adhered categories so that we can sum across
FD_New <- Final_Data %>%
  mutate(gro_1 = coalesce(gro_1, raft)) 
# This worked to fill in all N/A's in gro_1 column with 1 if there was a 1 in raft column

FD_New <- FD_New %>%
  mutate(gro_1 = coalesce(gro_1, adhered))
#Did this for gro column and adhered columns as well

#Now also fill in raft column for each adhered column
FD_New <- FD_New %>%
  mutate(raft = coalesce(raft, adhered))


#next  sum across growth columns to make the growth index value 
#This code is an if else to essentially set growth to 0 when the sample is arrested,
#also those replicates that transferred rafts that were MES and never grew we  keep at no growth
FD_New_t <- FD_New %>% rowwise() %>% 
  mutate(gro_index_t = case_when(
    is.na(arrested) ~ (sum(gro_1, raft, adhered, mot, ad_to_med, na.rm=TRUE)),
    arrested == 1 ~ 1))

#Then summarise/average growth across the replicates in treatment, sample group, and days
FD_Collapsed_t <- FD_New_t %>% group_by(treatment, samp_group, d_s_tran) %>%
  summarise(avg_grow= mean(gro_index_t),sd_grow = sd(gro_index_t))

#And again faceted for all treatment groups
FD_Collapsed_t$treatment_f = factor(FD_Collapsed_t$treatment, levels=c('A-Pos','F-Pos','P-Pos','P-Rinse',
                                                                       'P-dry short','P-dry long','P-NaCl low','P-NaCl high', "P-Neg"))
#####Plot for all treatments across time
FD_Collapsed_t %>%ggplot(aes(x=d_s_tran, y=avg_grow, 
                             color=as.factor(samp_group))) + 
  geom_point(position=position_jitter(w=0, h=0.1))+
  geom_line(position=position_jitter(w=0, h=0.1))+
  scale_x_continuous(breaks = seq(0, 30, by = 2))+
  geom_errorbar(aes(ymin=avg_grow-sd_grow, ymax=avg_grow+sd_grow), width=.2,
  position=position_dodge(0.05)) +
  facet_wrap(~treatment_f, ncol=3)

################################################################################
library(permute)
library(statmod)
library(plotrix)
#Proceeding with permutation test framework using package permute
#non parametric model; permutation test, separate tests for diff treatments.
#using two main metrics;
      #1. First comparing difference in mean GQI score between P-pos and treatments
      #2. Second comparing difference in viability between p-pos and treatments

#Analysis One: COMPARING MEANS
#Test 1: P_dry long compared to P-pos

target1 <- c("P-Pos","P-dry long")#Select treatments

meanDif <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-dry long"])
} #Create function for comparing means

Dry_long_test_mean <- FD_New_t %>% filter(d_s_tran < 15) %>% #Filter to dataset we want
  filter(samp_group < 15) %>%
  filter(treatment %in% target1) %>%
  select(treatment, gro_index_t, samp_group,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

#Compare means using two methods
diff_test<-Dry_long_test_mean %>% group_by(treatment) %>%
  summarise(grow = mean(grow), se = std.error(grow)) #Just a sumarise code to look at means of each group
diff_test
meanDif(Dry_long_test_mean$grow,Dry_long_test_mean$treatment)# TEST gives us a very similar answer to the with code below.

mtest_CTRL1 <- how(blocks=Dry_long_test_mean$samp_group) 

m_test1 <-shuffleSet(Dry_long_test_mean,nset=60000, control=mtest_CTRL1) 
numPerms(Dry_long_test_mean, control = mtest_CTRL1) #check number of possible perms with ctrl settings

#Run the permutation test
DDry_long_test_mean <- numeric(length = 60001)# this code creates the empty vector and number of perms
N <- nrow(Dry_long_test_mean)#Sets # of perms to length of empty vector
set.seed(42)# this makes it so randomization is the same each time, completely reproducible
for(i in seq_len(length(DDry_long_test_mean) - 1)) {
  perm <- shuffle(N, control=mtest_CTRL1) 
  DDry_long_test_mean[i] <- with(Dry_long_test_mean, meanDif(grow, treatment[perm]))
}# forloop to run perms, for each permutation, the difference in mean gro index for each treatment group is calculated and stored in the vector
DDry_long_test_mean[60001] <- with(Dry_long_test_mean, meanDif(grow, treatment))# fills in the last row in the vector with the mean diff from our observed data

#Plot for the perm test, showing created dist and our observed mean dif added
hist(DDry_long_test_mean, main = "",
     xlab = expression("Difference (P-Pos - P-dry long) in mean growth"))
rug(DDry_long_test_mean[60001], col = "red", lwd = 5) + 
  rug(QT_dry_long_mean, col = "blue", lwd =5)
#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot
#Also added 2.5% and 97.5% quartiles of the null distribution using the quantiles function below.

#Creating the p-value by counting the number of perms 
#where the differences in means were less than that seen in our observed data 
(Dbig1 <- sum(DDry_long_test_mean >= DDry_long_test_mean[60001]))

#Then dividing by total number of permutations #01/21/25 Updated using exact P-value formula for Monte Carlo based on Phipson and Smyth 2011
(Dbig1 +1) / (length(DDry_long_test_mean) + 1)

#using exact p-value methods from Phipson and Smyth #THIS IS USED FOR FINAL P-VALUES
permp(Dbig1,60001,11,12,method="exact",twosided = TRUE) 

#Getting 2.5% and 97.5% quantiles for histogram plots
QT_dry_long_mean <- quantile(DDry_long_test_mean, probs = c(0.025, 0.975), type=8)
QT_dry_long_mean

#Test Result 
#Diff in treatment means = 3.02; Pos mean 4.35, dry treatment mean 1.32
#Significant with p-value p < 0.001 (only 2 of 60001 showed)


#Test 2: P_dry long compared to P-pos
target2 <- c("P-Pos","P-Rinse")#Select treatments

meanDif2 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-Rinse"])
} #Create function for comparing means
Rinse_test_mean <- FD_New_t %>% filter(d_s_tran < 15) %>% #Filter to dataset we want
  filter(samp_group < 15) %>%
  filter(treatment %in% target2) %>%
  select(treatment, gro_index_t, samp_group,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

#Compare means using two methods
diff_test2<-Rinse_test_mean %>% group_by(treatment) %>%
  summarise(grow = mean(grow), se=std.error(grow)) #Just a summarize code to look at means of each group
diff_test2
meanDif2(Rinse_test_mean$grow,Rinse_test_mean$treatment)# TEST gives us a very similar answer to the with code below.

mtest_CTRL2 <- how(blocks=Rinse_test_mean$samp_group) 

m_test2 <-shuffleSet(Rinse_test_mean,nset=100, control=mtest_CTRL2) 
head(m_test2)
numPerms(Rinse_test_mean, control = mtest_CTRL2) #check number of possible perms with ctrl settings

#Run the permutation test
DRinse_test_mean <- numeric(length = 60001)# this code creates the empty vector and number of perms
N <- nrow(Rinse_test_mean)#Sets # of perms to length of empty vector
set.seed(123)# this makes it so randomization is the same each time, completely reproducible
for(i in seq_len(length(DRinse_test_mean) - 1)) {
  perm <- shuffle(N, control = mtest_CTRL2) 
  DRinse_test_mean[i] <- with(Rinse_test_mean, meanDif2(grow, treatment[perm]))
}# forloop to run perms, for each permutation, the difference in mean gro index for each treatment group is calculated and stored in the vector
DRinse_test_mean[60001] <- with(Rinse_test_mean, meanDif2(grow, treatment))# fills in the last row in the vector with the mean diff from our observed data

#Plot for the perm test, showing created dist and our observed mean dif added
hist(DRinse_test_mean, main = "",
     xlab = expression("Difference (P-Pos - P-Rinse) in mean growth"))
rug(DRinse_test_mean[60001], col = "red", lwd = 5) +
  rug(QT_rinse_mean, col = "blue", lwd =5)
#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

#Creating the p-value by counting the number of perms 
#where the differences in means were less than that seen in our observed data 
(Dbig2 <- sum(DRinse_test_mean >= DRinse_test_mean[60001]))

#Then dividing by total number of permutations #01/21/25 Updated using exact P-value formula for Monte Carlo based on Phipson and Smyth 2011
(Dbig2 +1) / (length(DRinse_test_mean) + 1)

permp(Dbig2,60001,11,12,method="exact",twosided = TRUE) #using exact p-value methods from Phipson and Smyth

#Getting 2.5% and 97.5% quantiles
QT_rinse_mean <- quantile(DRinse_test_mean, probs = c(0.025, 0.975), type=8)
QT_rinse_mean

#Results
#Diff in treatment means = 1.72; Pos mean 4.35, dry treatment mean 2.62
#Significant with p-value < 0.001 (4 of 60001)


#Test 3: P-NaCl high compared to P-pos
target3 <- c("P-Pos","P-NaCl high")#Select treatments

meanDif3 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-NaCl high"])
} #Create function for comparing means
Hnacl_test_mean <- FD_New_t %>% filter(d_s_tran < 15) %>% #Filter to dataset we want
  filter(samp_group < 15) %>%
  filter(treatment %in% target3) %>%
  select(treatment, gro_index_t, samp_group,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))
head(Hnacl_test_mean)
#Compare means using two methods
diff_test3<-Hnacl_test_mean %>% group_by(treatment) %>%
  summarise(grow = mean(grow), se=std.error(grow)) #Just a sumarise code to look at means of each group
diff_test3
meanDif3(Hnacl_test_mean$grow,Hnacl_test_mean$treatment)# TEST gives us a very similar answer to the with code below.

mtest_CTRL3 <- how(blocks=Hnacl_test_mean$samp_group) 

m_test3 <-shuffleSet(Hnacl_test_mean,nset=100, control=mtest_CTRL3) 
head(m_test3)
numPerms(Hnacl_test_mean, control = mtest_CTRL3)

#Run the permutation test
DHnacl_test_mean <- numeric(length = 60001)# this code creates the empty vector and number of perms
N <- nrow(Hnacl_test_mean)#Sets # of perms to length of empty vector
set.seed(123)# this makes it so randomization is the same each time, completely reproducible
for(i in seq_len(length(DHnacl_test_mean) - 1)) {
  perm <- shuffle(N, control=mtest_CTRL3) 
  DHnacl_test_mean[i] <- with(Hnacl_test_mean, meanDif3(grow, treatment[perm]))
}# forloop to run perms, for each permutation, the difference in mean gro index for each treatment group is calculated and stored in the vector
DHnacl_test_mean[60001] <- with(Hnacl_test_mean, meanDif3(grow, treatment))# fills in the last row in the vector with the mean diff from our observed data

#Plot for the perm test, showing created dist and our observed mean dif added
hist(DHnacl_test_mean, main = "",
     xlab = expression("Difference (P-Pos - P-NaCl high) in mean growth"))
rug(DHnacl_test_mean[60001], col = "red", lwd = 5) +
  rug(QT_Hnacl_mean, col ="blue", lwd=5)
#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

#Creating the p-value by counting the number of perms 
#where the differences in means were less than that seen in our observed data 
(Dbig3 <- sum(DHnacl_test_mean >= DHnacl_test_mean[60001]))

#Then dividing by total number of permutations 
(Dbig3 +1) / (length(DHnacl_test_mean) + 1)

permp(Dbig3,60001,11,12,method="exact",twosided = TRUE) #using exact p-value methods from Phipson and Smyth

#Getting 2.5% and 97.5% quantiles
QT_Hnacl_mean <- quantile(DHnacl_test_mean, probs = c(0.025, 0.975), type=8)
QT_Hnacl_mean

#Results
#Diff in treatment means = 0.79; Pos mean 4.35, dry treatment mean 3.56
#Marg significant p-value = .054


################################################################################
#Analysis Two: COMPARING VIABILITY

#Test 1: Long dry treatment pairwise comparison
target1 <- c("P-Pos","P-dry long")

FD_New_t$gro_index_t <- as.numeric(FD_New_t$gro_index_t)

Dry_long_test <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 15) %>%
  filter(treatment %in% target1) %>%
  select(treatment, gro_index_t, samp_group,rep) %>%
  group_by(treatment,samp_group, rep) %>%
  summarise(tot_obs = n(), grow = sum(gro_index_t > 1)) %>%
  mutate(grow_obs_ind = grow/tot_obs) #dplyr code to generate proportional viability metric

GOI_diff1 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-dry long"])
}    
GOI_diff1(Dry_long_test$grow_obs_ind,Dry_long_test$treatment)
Dry_long_test %>% group_by(treatment) %>%
  summarise(mean(grow_obs_ind))

rtest_CTRL1 <- how(blocks=Dry_long_test$samp_group) 
r_test1 <-shuffleSet(Dry_long_test,nset=100, control=rtest_CTRL1) 
head(r_test1)
numPerms(Dry_long_test, control = rtest_CTRL1)#check number of possible perms 

#Note not controlling for any structure. 
DDry_long_test <- numeric(length = 60001)# this code creates the empty vector and number of perms
N <- nrow(Dry_long_test)#Sets # of perms to length of empty vector
set.seed(13)# this makes it so randomization is the same each time, completely reproducible
for(i in seq_len(length(DDry_long_test) - 1)) {
  perm <- shuffle(N, control=rtest_CTRL1) 
  DDry_long_test[i] <- with(Dry_long_test, GOI_diff1(grow_obs_ind, treatment[perm]))
}# forloop to run perms, for each permutation, the difference in mean gro index for each treatment group is calulated and stored in the vector
DDry_long_test[60001] <- with(Dry_long_test, GOI_diff1(grow_obs_ind, treatment))# fills in the last row in the vector with the mean dif from our observed data

hist(DDry_long_test, main = "",
     xlab = expression("Difference (P-Pos - P-dry long) in mean viability"))
rug(DDry_long_test[60001], col = "red", lwd = 5) +
  rug(QT_drylong, col ="blue", lwd = 5)
  #plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

#Creating the "p-value" by counting the number of perms 
#where the differences in means were less than that seen in our observed data 
(Dbig5 <- sum(DDry_long_test>= DDry_long_test[60001]))

#Then dividing by total number of permutations #01/21/25 Updated using exact P-value formula for Monte Carlo based on Phipson and Smyth 2011
(Dbig5 +1)  / (length(DDry_long_test) +1)

permp(Dbig5,60001,11,12,method="exact",twosided = TRUE) #using exact p-value methods from Phipson and Smyth

#Code to generate 2.5 & 97.5 quantiles, using type 8 which is median unbiased regardless of distribution
QT_drylong <- quantile(DDry_long_test, probs = c(0.025, 0.975), type=8)
QT_drylong

#Results
#Diff in mean ratio = 56.9%; Pos ratio 100%, dry treatment ratio ~43.1%
#Significant with p-value = 0.0039.


#TEST 2 FOR RINSE
target2 <- c("P-Pos","P-Rinse")
FD_New_t$gro_index_t <- as.numeric(FD_New_t$gro_index_t)

Rinse_gro_thresh <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 15) %>%
  filter(treatment %in% target2) %>%
  select(treatment, gro_index_t, samp_group,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(tot_obs = n(), grow = sum(gro_index_t > 1)) %>%
  mutate(grow_obs_ind = grow/tot_obs)

ggplot(data=Rinse_gro_thresh, aes(x=as.factor(samp_group),y=grow_obs_ind)) +
  geom_col() +
  facet_wrap(~treatment, ncol=2)

GOI_diff2 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-Rinse"])
}    

GOI_diff2(Rinse_gro_thresh$grow_obs_ind,Rinse_gro_thresh$treatment)
Rinse_gro_thresh %>% group_by(treatment) %>%
  summarise(mean(grow_obs_ind))

rtest_CTRL2 <- how(blocks=Rinse_gro_thresh$samp_group) 
r_test2 <-shuffleSet(Rinse_gro_thresh,nset=100, control=rtest_CTRL2) 
head(r_test2)
numPerms(Rinse_gro_thresh, control = rtest_CTRL2)#check number of possible perms 

DRinse_gro_thresh <- numeric(length = 60001)# this code creates the empty vector and number of perms
N <- nrow(Rinse_gro_thresh)#Sets # of perms to length of empty vector
set.seed(13)# this makes it so randomization is the same each time, completely reproducible
for(i in seq_len(length(DRinse_gro_thresh) - 1)) {
  perm <- shuffle(N, control=rtest_CTRL2) 
  DRinse_gro_thresh[i] <- with(Rinse_gro_thresh, GOI_diff2(grow_obs_ind, treatment[perm]))
}# forloop to run perms, for each permutation, the difference in mean gro index for each treatment group is calulated and stored in the vector
DRinse_gro_thresh[60001] <- with(Rinse_gro_thresh, GOI_diff2(grow_obs_ind, treatment))# fills in the last row in the vector with the mean dif from our observed data

hist(DRinse_gro_thresh, main = "",
     xlab = expression("Difference (P-Pos - Rinse) in mean viability"))
rug(DRinse_gro_thresh[60001], col = "red", lwd = 5) +
  rug(QT_rinse_via, col ="blue", lwd = 5)
#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

#Creating the "p-value" by counting the number of perms 
#where the differences in means were less than that seen in our observed data 
(Dbig6 <- sum(DRinse_gro_thresh >= DRinse_gro_thresh[60001]))

#Then dividing by total number of permutations 
(Dbig6+1)  / (length(DRinse_gro_thresh)+1)

permp(Dbig6,60001,11,12,method="exact",twosided = TRUE) #using exact p-value methods from Phipson and Smyth

#Code to generate 2.5 & 97.5 quantiles, using type 8 which is median unbiased regardless of distribution
QT_rinse_via <- quantile(DRinse_gro_thresh, probs = c(0.025, 0.975), type=8)
QT_rinse_via

#Results
##Diff in mean ratio = 22.22%; Pos ratio 100%, rinse treatment ratio 77.8%
#Rinsing is significant reduction in observations over the growth threshold by 22%, p value = 0.001


#TEST 3: for high salt
target3 <- c("P-Pos","P-NaCl high")
FD_New_t$gro_index_t <- as.numeric(FD_New_t$gro_index_t)

NaClhigh_test <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 15) %>%
  filter(treatment %in% target3) %>%
  select(treatment, gro_index_t, samp_group,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(tot_obs = n(), grow = sum(gro_index_t > 1)) %>%
  mutate(grow_obs_ind = grow/tot_obs)

GOI_diff3 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-NaCl high"])
}    

GOI_diff3(NaClhigh_test$grow_obs_ind,NaClhigh_test$treatment)
NaClhigh_test %>% group_by(treatment) %>%
  summarise(mean(grow_obs_ind))

rtest_CTRL3 <- how(blocks=NaClhigh_test$samp_group) 
r_test3 <-shuffleSet(NaClhigh_test,nset=100, control=rtest_CTRL3) 
head(r_test3)

DNaClhigh_test <- numeric(length = 60001)# this code creates the empty vector and number of perms
N <- nrow(NaClhigh_test)#Sets # of perms to length of empty vector
set.seed(13)# this makes it so randomization is the same each time, completely reproducible
for(i in seq_len(length(DNaClhigh_test) - 1)) {
  perm <- shuffle(N, control=rtest_CTRL3) 
  DNaClhigh_test[i] <- with(NaClhigh_test, GOI_diff3(grow_obs_ind, treatment[perm]))
}# forloop to run perms, for each permutation, the difference in mean gro index for each treatment group is calulated and stored in the vector
DNaClhigh_test[60001] <- with(NaClhigh_test, GOI_diff3(grow_obs_ind, treatment))# fills in the last row in the vector with the mean dif from our observed data

hist(DNaClhigh_test, main = "",
     xlab = expression("Difference (P-Pos - High NaCl) in mean viability"))
rug(DNaClhigh_test[60001], col = "red", lwd = 5) +
  rug(QT_Hnacl_via, col ="blue", lwd = 3) #plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

#Creating the "p-value" by counting the number of perms 
#where the differences in means were less than that seen in our observed data 
(Dbig7 <- sum(DNaClhigh_test >= DNaClhigh_test[60001]))

#Then dividing by total number of permutations #01/21/25 Updated using exact P-value formula for Monte Carlo based on Phipson and Smyth 2011
(Dbig7 +1)  / (length(DNaClhigh_test) +1)

permp(Dbig7,60001,11,12,method="exact",twosided = TRUE) #using exact p-value methods from Phipson and Smyth

#Code to generate 2.5 & 97.5 quantiles, using type 8 which is median unbiased regardless of distribution
QT_Hnacl_via <- quantile(DNaClhigh_test, probs = c(0.025, 0.975), type=8)
QT_Hnacl_via

#Results
#Diff in mean ratio = 20.8%; Pos ratio 100%, rinse treatment ratio 79.1%
#Reduction not significant, p=.1

#TEST 4: short dry
target4 <- c("P-Pos","P-dry short")

dshort_test <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 15) %>%
  filter(treatment %in% target4) %>%
  select(treatment, gro_index_t, samp_group,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(tot_obs = n(), grow = sum(gro_index_t > 1)) %>%
  mutate(grow_obs_ind = grow/tot_obs)

GOI_diff4 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-dry short"])
}    

GOI_diff4(dshort_test$grow_obs_ind,dshort_test$treatment)
dshort_test %>% group_by(treatment) %>%
  summarise(mean(grow_obs_ind))

rtest_CTRL4 <- how(blocks=dshort_test$samp_group) 
r_test4 <-shuffleSet(dshort_test,nset=100, control=rtest_CTRL4) 
head(r_test4)

Ddshort_test <- numeric(length = 60001)# this code creates the empty vector and number of perms
N <- nrow(dshort_test)#Sets # of perms to length of empty vector
set.seed(13)# this makes it so randomization is the same each time, completely reproducible
for(i in seq_len(length(Ddshort_test) - 1)) {
  perm <- shuffle(N, control=rtest_CTRL4) #HOW DO WE KNOW SHUFFLE IS WORKING ON THE CORRECT WAY --> 01/16/2025 I checked and it is shuffling correctly
  Ddshort_test[i] <- with(dshort_test, GOI_diff4(grow_obs_ind, treatment[perm]))
}# forloop to run perms, for each permutation, the difference in mean gro index for each treatment group is calulated and stored in the vector
Ddshort_test[60001] <- with(dshort_test, GOI_diff4(grow_obs_ind, treatment))# fills in the last row in the vector with the mean dif from our observed data

hist(Ddshort_test, main = "",
     xlab = expression("Difference (P-Pos - short dry) in mean viability"))
rug(Ddshort_test[60001], col = "red", lwd = 5) +
  rug(QT_dryS_via, col= "blue", lwd = 3 )
  #plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

#Creating the "p-value" by counting the number of perms 
#where the differences in means were less than that seen in our observed data 
(Dbig8 <- sum(Ddshort_test >= Ddshort_test[60001]))

#Then dividing by total number of permutations 
(Dbig8 + 1)  / (length(Ddshort_test) +1)

permp(Dbig8,60001,11,12,method="exact",twosided = TRUE) #using exact p-value methods from Phipson and Smyth

#Code to generate 2.5 & 97.5 quantiles, using type 8 which is median unbiased regardless of distribution
QT_dryS_via <- quantile(Ddshort_test, probs = c(0.025, 0.975), type=8)
QT_dryS_via

#Results
##Diff in mean ratio = 12.5%; Pos ratio 100%, rinse treatment ratio 87.5%
#Reduction not significant, observations occurrs in 30% of data p = .149

######################################################################################################################################################
#POST-HOC Permutation tests to compare within treatment-day 
#data and within treatment data over time
# NOTE: Analysis uses imputation for missing positive control replicate in the day 1 treatment group
#We imputed the lowest growth value for the positive controls to be conservative.


#Analysis 1: Comparing treatments to positive controls in the day 1 treatment group

#TEST 1: Comparing p-pos to p-dry long treatment in day 1 treatment group
target1 <- c("P-Pos","P-dry long")

D1Drytest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group == 1) %>%
  filter(treatment %in% target1) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

#The following code imputes data in the positive control group
#adding the lowest observed growth value as a new row in the data   
imputed_row <- D1Drytest[2, ]
D1Drytest_I <- rbind(D1Drytest, imputed_row)
D1Drytest_I <- D1Drytest_I %>% arrange(desc(grow))
D1Drytest_I

meanDif1 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-dry long"])
}    

meanDif1(D1Drytest_I$grow,D1Drytest_I$treatment)
D1Drytest_I %>% group_by(treatment) %>%
  summarise(mean(grow))

ctrl_all <- how(complete = TRUE) 
allPerms(D1Drytest_I, control=ctrl_all)
numPerms(D1Drytest_I) #check number of possible perms with ctrl settings
Ph_shuff_test1 <-shuffleSet(D1Drytest_I, nset=120, control=ctrl_all) 
Ph_shuff_test1

#Note: Permutation tests are implemented with different code for the post hoc tests 
#in order to use all possible permutations of the data.

# Generate all permutations of the treatment labels
# allPerms returns a matrix where each row is a permutation of the treatment vector
permutations1 <- allPerms(D1Drytest_I$treatment)
permutations1

# Number of permutations
n_permutations1 <- nrow(permutations1)

# Initialize a vector to store the permuted mean differences
DD1Drytest <- numeric(length = n_permutations1 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic1 <- meanDif1(D1Drytest_I$grow, D1Drytest_I$treatment)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1Drytest) - 1)) {
  permuted_treatment1 <- permutations1[i, ]  # Get the i-th permutation of treatment
  DD1Drytest[i] <- with(D1Drytest_I, meanDif1(D1Drytest_I$grow, treatment[permuted_treatment1]))  # Store the permuted statistic
}

# Store the observed statistic in the last position
DD1Drytest[n_permutations1 + 1] <- observed_statistic1
DD1Drytest

hist(DD1Drytest, main = "",
     xlab = expression("Difference (P-Pos - P-dry long) in mean growth"))
rug(DD1Drytest[720], col = "red", lwd = 5) +
  rug(QT_drylong, col= "blue", lwd = 3 ) #plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigP_h1 <- sum(DD1Drytest >= DD1Drytest[720]))
P1<- (DbigP_h1)  / (length(DD1Drytest))
P1

p1 <- permp(DbigP_h1,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth
#Note: this p-value is calucalted without replacement 
p1

QT_drylong <- quantile(DD1Drytest, probs = c(0.025, 0.975), type=8)
QT_drylong

#Result
#Diff in mean GQI = 3; Pos GQI = 3.17, long dry treatment GQI = .167
#Reduction is significant p= 0.027


#TEST 2: Rinse to positive control in day 1
target2 <- c("P-Pos","P-Rinse")

D1Rintest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 2) %>%
  filter(treatment %in% target2) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

imputed_row <- D1Rintest[2, ]
D1Rintest <- rbind(D1Rintest, imputed_row)
D1Rintest <- D1Rintest %>% arrange(desc(treatment))
D1Rintest

meanDif2 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-Rinse"])
}    

meanDif2(D1Rintest$grow,D1Rintest$treatment)
D1Rintest %>% group_by(treatment) %>%
  summarise(mean(grow))

#All perms code for tests
permutations2 <- allPerms(D1Rintest)

# Number of permutations
n_permutations2 <- nrow(permutations2)

# Initialize a vector to store the permuted mean differences
DD1Rintest <- numeric(length = n_permutations2 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic2 <- meanDif2(D1Rintest$grow, D1Rintest$treatment)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1Rintest) - 1)) {
  permuted_treatment2 <- permutations2[i, ]  # Get the i-th permutation of treatment
  DD1Rintest[i] <- with(D1Rintest, meanDif2(D1Rintest$grow, treatment[permuted_treatment2]))  # Store the permuted statistic
}

# Store the observed statistic in the last position
DD1Rintest[n_permutations2 + 1] <- observed_statistic2
DD1Rintest

hist(DD1Rintest, main = "",xlim = c(-1, 1),
     xlab = expression("Difference (P-Pos - P-Rinse) in mean growth"))
rug(DD1Rintest[720], col = "red",lwd = 5) +
  rug(QT_rinse_d1, col= "blue", lwd = 3 ) #plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigP_h2 <- sum(DD1Rintest >= DD1Rintest[720]))
P2<- (DbigP_h2)  / (length(DD1Rintest))
P2

p2 <- permp(DbigP_h2,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth
p2

QT_rinse_d1<- quantile(DD1Rintest, probs = c(0.025, 0.975), type=8)
QT_rinse_d1

#Results
#Diff in mean GQI = .7222; Pos GQI = 3.17, long dry treatment GQI = 2.444
#Reduction is not significant p= 0.126


#TEST 3: High salt vs positive controls on day 1
target3 <- c("P-Pos","P-NaCl high")

D1NaclHtest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 2) %>%
  filter(treatment %in% target3) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

imputed_row <- D1NaclHtest[5, ]
D1NaclHtest <- rbind(D1NaclHtest, imputed_row)
D1NaclHtest <- D1NaclHtest %>% arrange(desc(grow))
D1NaclHtest

meanDif3 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-NaCl high"])
}    

meanDif3(D1NaclHtest$grow,D1NaclHtest$treatment)
D1NaclHtest %>% group_by(treatment) %>%
  summarise(mean(grow))

#All perms code for tests
permutations3 <- allPerms(D1NaclHtest)

# Number of permutations
n_permutations3 <- nrow(permutations3)

# Initialize a vector to store the permuted mean differences
DD1NaclHtest <- numeric(length = n_permutations3 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic3 <- meanDif3(D1NaclHtest$grow, D1NaclHtest$treatment)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1NaclHtest) - 1)) {
  permuted_treatment3 <- permutations3[i, ]  # Get the i-th permutation of treatment
  DD1NaclHtest[i] <- with(D1NaclHtest, meanDif3(D1NaclHtest$grow, treatment[permuted_treatment3]))  # Store the permuted statistic
}

# Store the observed statistic in the last position
DD1NaclHtest[n_permutations3 + 1] <- observed_statistic3

DD1NaclHtest

hist(DD1NaclHtest, main = "",
     xlab = expression("Difference (P-Pos - NaCl High) in mean growth"))
rug(DD1NaclHtest[720], col = "red",lwd = 5) +
  rug(QT_naclH_d1, col= "blue", lwd = 3 ) #plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigP_h3 <- sum(DD1NaclHtest >= DD1NaclHtest[720]))
P3<- (DbigP_h3)  / (length(DD1NaclHtest))
P3

p3 <- permp(DbigP_h3,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth
p3

QT_naclH_d1<- quantile(DD1NaclHtest, probs = c(0.025, 0.975), type=8)
QT_naclH_d1

#Results
#Diff in mean GQI = 2.77779; Pos GQI = 3.17, long dry treatment GQI = 0.389
#Reduction is significant p= 0.027

#TEST 4: Low salt compared to positive controls on Day 1
target4 <- c("P-Pos","P-NaCl low")

D1NaclLtest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 2) %>%
  filter(treatment %in% target4) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

imputed_row <- D1NaclLtest[5, ]
D1NaclLtest <- rbind(D1NaclLtest, imputed_row)
D1NaclLtest <- D1NaclLtest %>% arrange(desc(grow))
D1NaclLtest

meanDif4 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-NaCl low"])
}    

meanDif4(D1NaclLtest$grow,D1NaclLtest$treatment)
D1NaclLtest %>% group_by(treatment) %>%
  summarise(mean(grow))

#All perms code for tests
permutations4 <- allPerms(D1NaclLtest)

# Number of permutations
n_permutations4 <- nrow(permutations4)

# Initialize a vector to store the permuted mean differences
DD1NaclLtest <- numeric(length = n_permutations4 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic4 <- meanDif4(D1NaclLtest$grow, D1NaclLtest$treatment)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1NaclLtest) - 1)) {
  permuted_treatment4 <- permutations4[i, ]  # Get the i-th permutation of treatment
  DD1NaclLtest[i] <- with(D1NaclLtest, meanDif4(D1NaclLtest$grow, treatment[permuted_treatment4]))  # Store the permuted statistic
}

# Store the observed statistic in the last position
DD1NaclLtest[n_permutations4 + 1] <- observed_statistic4

DD1NaclLtest

hist(DD1NaclLtest, main = "",
     xlab = expression("Difference (P-Pos - NaCl Low) in mean growth"))
rug(DD1NaclLtest[720], col = "red", lwd = 4) +
  rug(QT_naclL_d1, col= "blue", lwd = 3 )
#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigP_h4 <- sum(DD1NaclLtest >= DD1NaclLtest[720]))
P4<- (DbigP_h4)  / (length(DD1NaclLtest))
P4

p4 <- permp(DbigP_h4,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth
p4

QT_naclL_d1<- quantile(DD1NaclLtest, probs = c(0.025, 0.975), type=8)
QT_naclL_d1

#Results
#Diff in mean GQI = 0.8889 ; Pos GQI = 3.17, long dry treatment GQI = 2.28
#Reduction is significant p= 0.027

#TEST 5: Short dry compared to controls in day 1
target5 <- c("P-Pos","P-dry short")

D1DryStest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 2) %>%
  filter(treatment %in% target5) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

imputed_row <- D1DryStest[2, ]
D1DryStest <- rbind(D1DryStest, imputed_row)
D1DryStest <- D1DryStest %>% arrange(treatment)
D1DryStest

meanDif5 <- function(x, grp) {
  mean(x[grp == "P-Pos" ]) - mean(x[grp == "P-dry short"])
}    

meanDif5(D1DryStest$grow,D1DryStest$treatment)
D1DryStest %>% group_by(treatment) %>%
  summarise(mean(grow))

#All perms code for tests
permutations5 <- allPerms(D1DryStest)

# Number of permutations
n_permutations5 <- nrow(permutations5)

# Initialize a vector to store the permuted mean differences
DD1DryStest <- numeric(length = n_permutations5 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic5 <- meanDif5(D1DryStest$grow, D1DryStest$treatment)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1DryStest) - 1)) {
  permuted_treatment4 <- permutations4[i, ]  # Get the i-th permutation of treatment
  DD1DryStest[i] <- with(D1DryStest, meanDif5(D1DryStest$grow, treatment[permuted_treatment4]))  # Store the permuted statistic
}

DD1DryStest[n_permutations5 + 1] <- observed_statistic5

DD1DryStest

hist(DD1DryStest, main = "",
     xlab = expression("Difference (P-Pos - P-dry short) in mean growth"))
rug(DD1DryStest[720], col = "red", lwd = 2) +
  rug(QT_drys_d1, col= "blue", lwd = 3 )#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigP_h5 <- sum(DD1DryStest >= DD1DryStest[720]))
P5<- (DbigP_h5)  / (length(DD1DryStest))
P5

p5 <- permp(DbigP_h5,720,3,3,method="exact",twosided = FALSE)
p5

QT_drys_d1<- quantile(DD1DryStest, probs = c(0.025, 0.975), type=8)
QT_drys_d1

#Results
#Diff in mean GQI = 0.1111; Pos GQI = 3.17, long dry treatment GQI = 3.06
#Reduction is not significant p= 0.47



#Post-Hoc Analysis 2: Within Dry-Long treatment x day 

#TEST 1: D1:D3 Long Dry tx comparison
Txtarg <- "P-dry long"

D1xD3LDryTest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 4) %>%
  filter(treatment %in% Txtarg) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

mLDryDif <- function(x, grp) {
  mean(x[grp == "1" ]) - mean(x[grp == "3"]) #samp_groups are factors remember
}    

mLDryDif(D1xD3LDryTest$grow, D1xD3LDryTest$samp_group)
D1xD3LDryTest %>% group_by(samp_group) %>%
  summarise(mean(grow))

#All perms code for tests
permutations6 <- allPerms(D1xD3LDryTest)

# Number of permutations
n_permutations6 <- nrow(permutations6)

# Initialize a vector to store the permuted mean differences
DD1xD3LDryTest <- numeric(length = n_permutations6 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic6 <- mLDryDif(D1xD3LDryTest$grow, D1xD3LDryTest$samp_group)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1xD3LDryTest) - 1)) {
  permuted_treatment6 <- permutations6[i, ]  # Get the i-th permutation of treatment
  DD1xD3LDryTest[i] <- with(D1xD3LDryTest, mLDryDif(D1xD3LDryTest$grow, samp_group[permuted_treatment6]))  # Store the permuted statistic
}

DD1xD3LDryTest[n_permutations6 + 1] <- observed_statistic6

DD1xD3LDryTest

hist(DD1xD3LDryTest, main = "",
     xlab = expression("Difference (Long dry day1 - Long dry day3) in mean growth"))
rug(DD1xD3LDryTest[720], col = "red", lwd = 5) +
  rug(QT_D1XD3_via, col ="blue",lwd = 3)#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigPTx1 <- sum(DD1xD3LDryTest <= DD1xD3LDryTest[720]))
TxP1_3<- (DbigPTx1)  / (length(DD1xD3LDryTest))
TxP1_3

permp(DbigPTx1,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth

QT_D1XD3_via <- quantile(DD1xD3LDryTest, probs = c(0.025, 0.975), type=8)
QT_D1XD3_via

#Results
#Diff in mean GQI = 2.83333; Day1 GQI = 0.167, Day3 GQI = 3
#Reduction is significant p= 0.027

#TEST 2: D1:D7 long dry Tx compare
D1xD7LDryTest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group =="1"|samp_group =="7") %>%
  filter(treatment %in% Txtarg) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

mLDryDif1X7 <- function(x, grp) {
  mean(x[grp == "1" ]) - mean(x[grp == "7"]) #samp_groups are factors remember
}    

mLDryDif1X7 (D1xD7LDryTest$grow, D1xD7LDryTest$samp_group)
D1xD7LDryTest %>% group_by(samp_group) %>%
  summarise(mean(grow))

#All perms code for tests
permutations7 <- allPerms(D1xD7LDryTest)

# Number of permutations
n_permutations7 <- nrow(permutations7)

# Initialize a vector to store the permuted mean differences
DD1xD7LDryTest <- numeric(length = n_permutations7 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic7 <- mLDryDif1X7(D1xD7LDryTest$grow, D1xD7LDryTest$samp_group)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1xD7LDryTest) - 1)) {
  permuted_treatment7 <- permutations7[i, ]  # Get the i-th permutation of treatment
  DD1xD7LDryTest[i] <- with(D1xD7LDryTest, mLDryDif1X7(D1xD7LDryTest$grow, samp_group[permuted_treatment7]))  # Store the permuted statistic
}

DD1xD7LDryTest[n_permutations7 + 1] <- observed_statistic7

DD1xD7LDryTest


hist(DD1xD7LDryTest, main = "",
     xlab = expression("Difference (Long dry day1 - Long dry day7) in mean growth"))
rug(DD1xD7LDryTest[720], col = "red", lwd = 5) +
  rug(QT_D1XD7_via, col ="blue", lwd=3)#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigPTx2 <- sum(DD1xD7LDryTest <= DD1xD7LDryTest[720]))
TxP1_7<- (DbigPTx2)  / (length(DD1xD7LDryTest))
TxP1_7

permp(DbigPTx2,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth

QT_D1XD7_via <- quantile(DD1xD7LDryTest, probs = c(0.025, 0.975), type=8)
QT_D1XD7_via

#Results
#Diff in mean GQI = 1.2777; Day1 GQI = 0.167, Day3 GQI = 1.44
#Reduction is not significant different p= 0.47

#TEST 3: D1:D14 long dry Tx compare
D1xD14LDryTest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group =="1"|samp_group =="14") %>%
  filter(treatment %in% Txtarg) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

mLDryDif1X14 <- function(x, grp) {
  mean(x[grp == "1" ]) - mean(x[grp == "14"]) #samp_groups are factors remember
}    

mLDryDif1X14(D1xD14LDryTest$grow, D1xD14LDryTest$samp_group)
D1xD14LDryTest %>% group_by(samp_group) %>%
  summarise(mean(grow))

#All perms code for tests
permutations8 <- allPerms(D1xD14LDryTest)

# Number of permutations
n_permutations8 <- nrow(permutations8)

# Initialize a vector to store the permuted mean differences
DD1xD14LDryTest <- numeric(length = n_permutations8 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic8 <- mLDryDif1X14(D1xD14LDryTest$grow, D1xD14LDryTest$samp_group)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1xD14LDryTest) - 1)) {
  permuted_treatment8 <- permutations8[i, ]  # Get the i-th permutation of treatment
  DD1xD14LDryTest[i] <- with(D1xD14LDryTest, mLDryDif1X14(D1xD14LDryTest$grow, samp_group[permuted_treatment8]))  # Store the permuted statistic
}

DD1xD14LDryTest[n_permutations8 + 1] <- observed_statistic8

DD1xD14LDryTest

hist(DD1xD14LDryTest, main = "",
     xlab = expression("Difference (Long dry day1 - Long dry day14) in mean growth"))
rug(DD1xD14LDryTest[720], col = "red", lwd = 5) +
  rug(QT_D1XD14_via, col ="blue", lwd=3) #plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigPTx3 <- sum(DD1xD14LDryTest <= DD1xD14LDryTest[720]))
TxP1_14<- (DbigPTx3)  / (length(DD1xD14LDryTest))
TxP1_14

permp(DbigPTx3,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth

QT_D1XD14_via <- quantile(DD1xD14LDryTest, probs = c(0.025, 0.975), type=8)
QT_D1XD14_via

#Results
#Diff in mean GQI = 0.5 ; Day1 GQI = 0.167, Day3 GQI = 0.667
#Reduction is not significant different p= 0.47

###########
#Post-Hoc Analysis 3: Within high salt treatment x day 
#Use subset to only compare D1 to other days

#TEST 1: D1:D3 high salt tx comparison
TxHtarg <- "P-NaCl high"

D1xD3HNaCLTest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 4) %>%
  filter(treatment %in% TxHtarg) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

mHNaCLDif <- function(x, grp) {
  mean(x[grp == "1" ]) - mean(x[grp == "3"]) #samp_groups are factors remember
}    

mHNaCLDif(D1xD3HNaCLTest$grow, D1xD3HNaCLTest$samp_group)
D1xD3HNaCLTest %>% group_by(samp_group) %>%
  summarise(mean(grow))

#All perms code for tests
permutations9 <- allPerms(D1xD3HNaCLTest)

# Number of permutations
n_permutations9 <- nrow(permutations9)

# Initialize a vector to store the permuted mean differences
DD1xD3HNaCLTest <- numeric(length = n_permutations9 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic9 <- mHNaCLDif(D1xD3HNaCLTest$grow, D1xD3HNaCLTest$samp_group)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1xD3HNaCLTest) - 1)) {
  permuted_treatment9 <- permutations9[i, ]  # Get the i-th permutation of treatment
  DD1xD3HNaCLTest[i] <- with(D1xD3HNaCLTest, mHNaCLDif(D1xD3HNaCLTest$grow, samp_group[permuted_treatment9]))  # Store the permuted statistic
}

DD1xD3HNaCLTest[n_permutations9 + 1] <- observed_statistic9

DD1xD3HNaCLTest

hist(DD1xD3HNaCLTest, main = "",
     xlab = expression("Difference (high NaCl day1 - high NaCl day3) in mean growth"))
rug(DD1xD3HNaCLTest[720], col = "red", lwd = 5)+
  rug(QT_D1XD3_via, col= "blue", lwd= 3)#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigPTx1 <- sum(DD1xD3HNaCLTest <= DD1xD3HNaCLTest[720]))
TxP1_3<- (DbigPTx1)  / (length(DD1xD3HNaCLTest))
TxP1_3

permp(DbigPTx1,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth

# D1 Mean GQI was 2.833 lower than D3 mean GQI, significant lower p=.049
QT_D1XD3_via <- quantile(DD1xD3HNaCLTest, probs = c(0.025, 0.975), type=8)
QT_D1XD3_via

#Results
#Diff in mean GQI = 3.7778 ; Day1 GQI = 0.389, Day3 GQI = 4.17
#Reduction is significantly different p= 0.027


#TEST 2: D1:D7 high salt Tx compare
D1xD7HNaCLTest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group =="1"|samp_group =="7") %>%
  filter(treatment %in% TxHtarg) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

mHNaCLDif <- function(x, grp) {
  mean(x[grp == "1" ]) - mean(x[grp == "7"]) #samp_groups are factors remember
}    

mHNaCLDif(D1xD7HNaCLTest$grow, D1xD7HNaCLTest$samp_group)
D1xD7HNaCLTest %>% group_by(samp_group) %>%
  summarise(mean(grow))

#All perms code for tests
permutations10 <- allPerms(D1xD7HNaCLTest)

# Number of permutations
n_permutations10 <- nrow(permutations10)

# Initialize a vector to store the permuted mean differences
DD1xD7HNaCLTest <- numeric(length = n_permutations10 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic10 <- mHNaCLDif(D1xD7HNaCLTest$grow, D1xD7HNaCLTest$samp_group)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1xD7HNaCLTest) - 1)) {
  permuted_treatment10 <- permutations10[i, ]  # Get the i-th permutation of treatment
  DD1xD7HNaCLTest[i] <- with(D1xD7HNaCLTest, mHNaCLDif(D1xD7HNaCLTest$grow, samp_group[permuted_treatment10]))  # Store the permuted statistic
}

DD1xD7HNaCLTest[n_permutations10 + 1] <- observed_statistic10

DD1xD7HNaCLTest

hist(DD1xD7HNaCLTest, main = "",
     xlab = expression("Difference (high NaCl day1 - high NaCl day7) in mean growth"))
rug(DD1xD7HNaCLTest[720], col = "red", lwd = 5) +
  rug(QT_D1XD7_via, col ="blue", lwd=3)#plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigPTx2 <- sum(DD1xD7HNaCLTest <= DD1xD7HNaCLTest[720]))
TxP1_7<- (DbigPTx2)  / (length(DD1xD7HNaCLTest))
TxP1_7

permp(DbigPTx2,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth

QT_D1XD7_via <- quantile(DD1xD7HNaCLTest, probs = c(0.025, 0.975), type=8)
QT_D1XD7_via

#Results
#Diff in mean GQI = 3.77778 ; Day1 GQI = 0.389, Day3 GQI = 4.17
#Reduction is significantly different p= 0.027

#TEST 3: D1:D14 high salt Tx compare
D1xD14HNaCLTest <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group =="1"|samp_group =="14") %>%
  filter(treatment %in% TxHtarg) %>%
  select(treatment, gro_index_t, samp_group,d_s_tran,rep) %>%
  group_by(samp_group,rep) %>%
  summarise(grow = mean(gro_index_t))

mLDryDif1X14 <- function(x, grp) {
  mean(x[grp == "1" ]) - mean(x[grp == "14"]) #samp_groups are factors remember
}    

mLDryDif1X14(D1xD14HNaCLTest$grow, D1xD14HNaCLTest$samp_group)
D1xD14HNaCLTest %>% group_by(samp_group) %>%
  summarise(mean(grow))

#All perms code for tests
permutations11 <- allPerms(D1xD14HNaCLTest)

# Number of permutations
n_permutations11 <- nrow(permutations11)

# Initialize a vector to store the permuted mean differences
DD1xD14HNaCLTest <- numeric(length = n_permutations11 + 1)  # Add 1 for the observed statistic

# Calculate the observed mean difference
observed_statistic11 <- mLDryDif1X14(D1xD14HNaCLTest$grow, D1xD14HNaCLTest$samp_group)

# Loop over all permutations and calculate the mean difference for each
for (i in seq_len(length(DD1xD14HNaCLTest) - 1)) {
  permuted_treatment11 <- permutations11[i, ]  # Get the i-th permutation of treatment
  DD1xD14HNaCLTest[i] <- with(D1xD14HNaCLTest, mLDryDif1X14(D1xD14HNaCLTest$grow, samp_group[permuted_treatment11]))  # Store the permuted statistic
}

DD1xD14HNaCLTest[n_permutations11 + 1] <- observed_statistic11

DD1xD14HNaCLTest

hist(DD1xD14HNaCLTest, main = "",
     xlab = expression("Difference (high NaCl day1 - high NaCl day14) in mean growth"))
rug(DD1xD14HNaCLTest[720], col = "red", lwd = 5) +
  rug(QT_D1XD14_via, col="blue", lwd=3) #plots a histogram of the the created dist of mean differences, 
#and then plots our observed difference in a red rug plot

(DbigPTx3 <- sum(DD1xD14HNaCLTest <= DD1xD14HNaCLTest[720]))
TxP1_14<- (DbigPTx3)  / (length(DD1xD14HNaCLTest))
TxP1_14

permp(DbigPTx3,720,3,3,method="exact",twosided = FALSE) #using exact p-value methods from Phipson and Smyth

QT_D1XD14_via <- quantile(DD1xD14HNaCLTest, probs = c(0.025, 0.975), type=8)
QT_D1XD14_via

#Results
#Diff in mean GQI = 4.444 ; Day1 GQI = 0.389, Day3 GQI = 4.83
#Reduction is significantly different p= 0.027

################################################################################
#Figure Code
library(RColorBrewer)
library(grid)

#Organizing data for figures
FD_New_t_fig <- FD_New_t %>% mutate(treatment_name = 
                                      case_when(treatment == "P-Pos"
                                                ~ "Pelt positive",
                                                treatment == "P-NaCl high" ~ "High NaCl",
                                                treatment == "P-NaCl low" ~ "Low NaCl",
                                                treatment == "P-dry long" ~ "Long dry",
                                                treatment == "P-dry short" ~ "Short dry",
                                                treatment == "P-Rinse" ~ "Rinse",
                                                treatment == "A-Pos" ~ "Agar positive",
                                                treatment == "F-Pos" ~ "Cotton positive",
                                                treatment == "P-Neg" ~ "Pelt negative")) #Rename treatment groups for plot

FD_New_t_fig$treatment_name = factor(FD_New_t_fig$treatment_name, levels=c('Agar positive','Cotton positive','Pelt positive','Rinse',
                                                                           'Short dry','Long dry','Low NaCl','High NaCl', "Pelt negative"))
FD_New_t_fig  

#Figure 3 panel B
Diffplot2<- FD_New_t_fig %>% filter(d_s_tran < 15) %>% 
  filter(samp_group <= 14) %>% group_by(treatment,samp_group, rep) %>% 
  summarise(grow= mean(gro_index_t), se=std.error(gro_index_t)) #Gives us values

Diffplot3<- Diffplot2 %>% group_by(treatment, samp_group) %>% 
  summarise(growM= mean(grow), seM=std.error(grow)) #Gives us values

diff_from_control <- Diffplot3 %>% group_by(samp_group) %>%
  mutate(control_value = growM[treatment == "P-Pos"],  # Get the control value for each day
          value_diff = growM - control_value, # Calculate the difference with the control
         u_lim = (growM+seM) - control_value,
         l_lim = (growM-seM) - control_value,
         )%>%                                   
  ungroup() %>%
  select(-control_value)#Code to generate the diff column, partially chatGPT. 
#Could also add and subtract the SD to grow values to create error bars for this plot
#Add full treatment means
DiffplotTreatmentTot<- Diffplot2 %>%  group_by(treatment) %>% 
  summarise(growTot= mean(grow), seTot=std.error(grow)) #Gives us values

diff_from_control_Tot <- DiffplotTreatmentTot %>% 
  mutate(control_value = growTot[treatment == "P-Pos"],  # Get the control value for each day
         value_diff = growTot - control_value, # Calculate the difference with the control
         u_lim = (growTot+seTot) - control_value,
         l_lim = (growTot-seTot) - control_value,
  )%>%                                   
  ungroup() %>%
  select(-control_value)

diff_from_control_TotP <- diff_from_control_Tot %>% mutate(treatment_name = 
                               case_when(treatment == "P-Pos"
                                         ~ "Pelt positive",
                                         treatment == "P-NaCl high" ~ "High NaCl",
                                         treatment == "P-NaCl low" ~ "Low NaCl",
                                         treatment == "P-dry long" ~ "Long dry",
                                         treatment == "P-dry short" ~ "Short dry",
                                         treatment == "P-Rinse" ~ "Rinse",
                                         treatment == "A-Pos" ~ "Agar positive",
                                         treatment == "F-Pos" ~ "Cotton positive",
                                         treatment == "P-Neg" ~ "Pelt negative"))%>% 
  filter(treatment != "P-Neg") %>%
  filter(treatment != "P-Pos") %>%
  filter(treatment != "A-Pos") %>%
  filter(treatment != "F-Pos")

diff_from_control$samp_group <- as.factor(diff_from_control$samp_group)
diff_from_controlP <- diff_from_control %>%  
mutate(treatment_name = 
         case_when(treatment == "P-Pos"
                   ~ "Pelt positive",
                   treatment == "P-NaCl high" ~ "High NaCl",
                   treatment == "P-NaCl low" ~ "Low NaCl",
                   treatment == "P-dry long" ~ "Long dry",
                   treatment == "P-dry short" ~ "Short dry",
                   treatment == "P-Rinse" ~ "Rinse",
                   treatment == "A-Pos" ~ "Agar positive",
                   treatment == "F-Pos" ~ "Cotton positive",
                   treatment == "P-Neg" ~ "Pelt negative"))%>%
  filter(treatment != "P-Neg") %>%
  filter(treatment != "P-Pos") %>%
  filter(treatment != "A-Pos") %>%
  filter(treatment != "F-Pos")

diff_from_controlP$samp_group <- as.factor(diff_from_controlP$samp_group)
diff_from_controlP<- diff_from_controlP %>% add_row(treatment_name="Rinse", 
                                    samp_group = "treatment\nmean") 


diff_from_controlP$samp_group <- as.factor(diff_from_controlP$samp_group)
diff_from_controlP$samp_group<- factor(diff_from_controlP$samp_group,levels=
                                  c('1','3',
                                    '7','14',
                                    'treatment\nmean'))

diff_from_controlP$treatment_name = factor(diff_from_controlP$treatment_name,
                                           levels=c('Rinse','Short dry','Long dry',
                                                    'Low NaCl','High NaCl','treatment\nmean'))
#Code above sets treatment names for plot and removes negative and positive controls

#this code sets jitters of the data so that I can jitter the SE bars and points the same way

diff_from_control_TotP <- diff_from_control_TotP %>%
  mutate(jitter_x = 5 + rnorm(n(), 0.1, 0.1)) 

diff_from_control_plot <- diff_from_controlP %>% arrange(treatment_name, samp_group,) %>%
  ggplot(aes(x=samp_group, y=value_diff, 
          color=as.factor(treatment_name), group=as.factor(treatment_name))) + 
  geom_point(position = position_dodge(width = 0.3)
             ,size=3.5) +
  geom_line(position = position_dodge(width = 0.3), 
            size=.75, alpha =0.35)+
  geom_linerange(aes(x = as.factor(samp_group), 
                     ymin=l_lim, ymax =u_lim,
                     color=as.factor(treatment_name)),
                 position = position_dodge(width = 0.3),
                 size=.75,) +
  labs(x="Treatment day", y=expression("Difference in growth from control"),
       color="Treatment") +
  labs(shape = element_blank()) +
  theme_bw() +
  geom_abline(slope=0, intercept=0, linetype="dashed", color = "black") +
  geom_point(data=diff_from_control_TotP, 
             aes(x =jitter_x, y=value_diff,
                 color=treatment_name),size=3.5,
             inherit.aes = FALSE) +
  geom_linerange(data=diff_from_control_TotP,
                 aes(x = jitter_x,
                     ymin=l_lim, ymax =u_lim,
                 color=treatment_name), size=.75, 
                 inherit.aes = FALSE) + 
  geom_vline(xintercept = 4.4, color="gray50",size=.75)+
  scale_color_brewer(palette = "Set1") 

asterisk_annotation1 <- textGrob("*", gp = gpar(fontsize = 19, col = "black"))
asterisk_annotation2 <- textGrob("*", gp = gpar(fontsize = 19, col = "black"))

 diff_f_con_FIN <- diff_from_control_plot + 
  annotation_custom(asterisk_annotation1, xmin = 5.34, xmax = 5.34, 
                    ymin = -3, ymax = -3) +
  annotation_custom(asterisk_annotation2, xmin = 5.13, xmax = 5.13, 
                    ymin = -1.7, ymax = -1.7)

diff_grow_pan<- diff_f_con_FIN + theme(axis.text=element_text(size=12),
                        axis.title=element_text(size=14),
                        legend.title=element_text(size=14), 
                        legend.text=element_text(size=12))

diff_grow_pan #Final panel 

#Save code
ggsave(filename="diff_f_con_FIN_plot.png", plot= diff_f_con_FIN , 
       device="png",
       path="C:/Users/jabur/OneDrive - Washington State University (email.wsu.edu)/Documents/Fomite Study Materials/WYGF-CPW Experiment/April 29 Experiment/R_Data/Pelt_Bd",
       height=4, width=7, units="in", dpi=800)


#Figure 3 panel A: Diff from control viability
head(FD_New_t)

ratio_plot <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 15) %>%
  select(treatment, gro_index_t, samp_group, rep) %>%
  group_by(treatment,samp_group,rep) %>%
  summarise(tot_obs = n(), grow = sum(gro_index_t > 1)) %>%
  mutate(gro_ratio = grow/tot_obs)
head(ratio_plot)

ratio_plot_byrep <-ratio_plot %>%
  group_by(treatment,samp_group) %>%
  summarize(ratio_mean = mean(gro_ratio), ratio_se = std.error(gro_ratio))

ratio_diff_cont <- ratio_plot_byrep %>% group_by(samp_group) %>%
  mutate(control_value = ratio_mean[treatment == "P-Pos"],  # Get the control value for each day
         value_diff = ratio_mean - control_value,
         u_lim = (ratio_mean + ratio_se) - control_value,
         l_lim = (ratio_mean - ratio_se) - control_value# Calculate the difference with the control
)%>%                                   
  ungroup() %>%
  select(-control_value)

##This code for making treatment group totals across all days
ratio_diff_TOT<- ratio_plot %>%
  group_by(treatment) %>%
  summarise(TOT = sum(tot_obs), TOT_Grow = sum(grow)) %>%
  mutate(TOT_Ratio = TOT_Grow/TOT) %>%
  mutate(TOT_Diff = TOT_Ratio -1)

ratio_diff_TOT <-ratio_plot %>%
  group_by(treatment) %>%
  summarize(TOT_mean = mean(gro_ratio), TOT_se = std.error(gro_ratio))

diff_from_control_Tot <- ratio_diff_TOT %>% 
  mutate(control_value = TOT_mean[treatment == "P-Pos"],  # Get the control value for each day
         value_diff = TOT_mean - control_value, # Calculate the difference with the control
         u_lim = (TOT_mean+TOT_se) - control_value,
         l_lim = (TOT_mean-TOT_se) - control_value,
  )%>%                                   
  ungroup() %>%
  select(-control_value)

##Code above sets treatment names for plot and removes negative and positive controls
ratio_diff_contP <- ratio_diff_cont %>% mutate(treatment_name = 
                                                     case_when(treatment == "P-Pos"
                                                               ~ "Pelt positive",
                                                               treatment == "P-NaCl high" ~ "High NaCl",
                                                               treatment == "P-NaCl low" ~ "Low NaCl",
                                                               treatment == "P-dry long" ~ "Long dry",
                                                               treatment == "P-dry short" ~ "Short dry",
                                                               treatment == "P-Rinse" ~ "Rinse",
                                                               treatment == "A-Pos" ~ "Agar positive",
                                                               treatment == "F-Pos" ~ "Cotton positive",
                                                               treatment == "P-Neg" ~ "Pelt negative"))%>% 
  filter(treatment != "P-Neg") %>%
  filter(treatment != "P-Pos") %>%
  filter(treatment != "A-Pos") %>%
  filter(treatment != "F-Pos")

ratio_diff_contP$samp_group <- as.factor(ratio_diff_contP$samp_group) 

ratio_diff_contP<- ratio_diff_contP %>% add_row(treatment_name="Rinse", 
                                                samp_group = "treatment\nmean")

ratio_diff_contP$treatment_name = factor(ratio_diff_contP$treatment_name,
                                           levels=c('Rinse','Short dry','Long dry','Low NaCl','High NaCl'))

ratio_diff_contP$samp_group<- factor(ratio_diff_contP$samp_group,levels=
                                      c('1','3',
                                        '7','14',
                                        'treatment\nmean'))

#Now for the treatment total data
ratio_diff_TOTP <- diff_from_control_Tot %>% mutate(treatment_name = 
                                                 case_when(treatment == "P-Pos"
                                                           ~ "Pelt positive",
                                                           treatment == "P-NaCl high" ~ "High NaCl",
                                                           treatment == "P-NaCl low" ~ "Low NaCl",
                                                           treatment == "P-dry long" ~ "Long dry",
                                                           treatment == "P-dry short" ~ "Short dry",
                                                           treatment == "P-Rinse" ~ "Rinse",
                                                           treatment == "A-Pos" ~ "Agar positive",
                                                           treatment == "F-Pos" ~ "Cotton positive",
                                                           treatment == "P-Neg" ~ "Pelt negative"))%>% 
  filter(treatment != "P-Neg") %>%
  filter(treatment != "P-Pos") %>%
  filter(treatment != "A-Pos") %>%
  filter(treatment != "F-Pos")

ratio_diff_TOTP$treatment_name = factor(ratio_diff_TOTP$treatment_name,
                                         levels=c('Rinse','Short dry','Long dry','Low NaCl','High NaCl'))

ratio_diff_TOTP  <- ratio_diff_TOTP  %>%
  mutate(jitter_x = 5 + rnorm(n(), 0.1, 0.1)) 

ratio_diff_contP <- ratio_diff_contP  %>%
  mutate(jitter_x = samp_group + rnorm(n(), 0.1, 0.1)) 

###PLOT CODE
ratio_diff_cont_plot <- ratio_diff_contP %>%
  ggplot(aes(x=as.factor(samp_group), y=value_diff, 
             color=as.factor(treatment_name), group=as.factor(treatment_name))) + 
  geom_point(position = position_dodge(width = 0.3)
             ,size=3.5) +
  geom_line(position = position_dodge(width = 0.3),size=.75, alpha =0.35)+
  geom_linerange(aes(x = as.factor(samp_group), 
                     ymin=l_lim, ymax =u_lim,
                     color=as.factor(treatment_name)),
                 position = position_dodge(width = 0.3),
                 size=.75,) +
  labs(x="Treatment day", y=expression("Difference in viability from control"),
       color="Treatment") +
  labs(shape = element_blank()) +
  theme_bw() +
  geom_abline(slope=0, intercept=0, linetype="dashed", color = "black") +
  geom_point(data=ratio_diff_TOTP, 
             aes(x =jitter_x, y=value_diff,
                 color=treatment_name),size=3.5,
             inherit.aes = FALSE)+ 
  geom_linerange(data=ratio_diff_TOTP,
                 aes(x = jitter_x,
                     ymin=l_lim, ymax =u_lim,
                     color=treatment_name), size=.75, 
                 inherit.aes = FALSE) +
  geom_vline(xintercept = 4.4, color="gray50",size=.75)+
  scale_color_brewer(palette = "Set1") 

asterisk_annotation3 <- textGrob("*", gp = gpar(fontsize = 19, col = "black"))
asterisk_annotation4 <- textGrob("*", gp = gpar(fontsize = 19, col = "black"))

ratio_diff_cont_plot_FIN  <- ratio_diff_cont_plot + 
  annotation_custom(asterisk_annotation3, xmin = 5.15, xmax = 5.15, 
                    ymin = -.565, ymax = -.565) +
  annotation_custom(asterisk_annotation4, xmin = 5.23, xmax = 5.23, 
                    ymin = -.22, ymax = -.22)

diff_ratio_pan<- ratio_diff_cont_plot_FIN + theme(axis.text=element_text(size=12),
                       axis.title=element_text(size=14),
                       legend.title=element_text(size=14), 
                       legend.text=element_text(size=12))
diff_ratio_pan #Final panel code

ggsave(filename="ratio_diff_cont_plot.png", plot=ratio_diff_cont_plot_FIN  , 
       device="png",
       path="C:/Users/jabur/OneDrive - Washington State University (email.wsu.edu)/Documents/Fomite Study Materials/WYGF-CPW Experiment/April 29 Experiment/R_Data/Pelt_Bd",
       height=4, width=7, units="in", dpi=800)

#Combine Figure 3 Panels into one plot
library(gridExtra)
library(cowplot)
plot1<- diff_ratio_pan + theme(axis.title.x = element_blank(),
                               legend.position="none",plot.margin = margin(t = 20)) 
plot2 <- diff_grow_pan + theme(legend.position = "bottom", 
                               legend.title = element_blank())


comb_plot<- plot_grid(plot1, plot2, ncol = 1, rel_heights = c(1, 1.2), align="v", axis="tb") 
comb_plot

ggsave(filename="comb_plot.png", plot=comb_plot  , 
       device="png",
       path="C:/Users/jabur/OneDrive - Washington State University (email.wsu.edu)/Documents/Fomite Study Materials/WYGF-CPW Experiment/April 29 Experiment/R_Data/Pelt_Bd",
       height=7.5, width=7, units="in", dpi=800)


#Supporting Information and exploratory plots
#Figure S19
A_GQI_ALLOBS<- FD_New_t_fig %>%ggplot(aes(x=d_s_tran, y=gro_index_t, 
                                          color=as.factor(samp_group))) + 
  stat_summary(position=position_dodge(width=1), size=.5)+
  stat_summary(fun.y=mean,position=position_dodge(width=2), geom="line", size=.5)+
  scale_x_continuous(breaks = seq(0, 30, by = 2))+
  labs(x="Observation Period (Days)", y="Average Growth Quality Index",
       color="Treatment day\n Group") +
  labs(shape = element_blank()) +
  theme_bw() +
  geom_abline(slope=0, intercept=1, linetype="dashed", color = "red") +
  geom_vline(xintercept=14, linetype="dashed", color = "black") +
  scale_y_continuous(breaks=c(0,1,2,3,4,5)) +
  coord_cartesian(ylim = c(0, 5))+
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30))+
  facet_wrap(~treatment_name, ncol=3)
A_GQI_ALLOBS

#plt save code
ggsave(filename="A_GQI_ALLOBS.png", plot=A_GQI_ALLOBS, device="png",
       path="C:/Users/jabur/OneDrive - Washington State University (email.wsu.edu)/Documents/Fomite Study Materials/WYGF-CPW Experiment/April 29 Experiment/R_Data/Pelt_Bd",
       height=7, width=10, units="in", dpi=500)

#Note: this figure does not appear in the text or supporting information
ggplot(FD_New_t, aes(x = gro_index_t)) +
  geom_histogram()+
  facet_wrap(~treatment, ncol=4)

Tot_grothresh <- FD_New_t %>% filter(d_s_tran < 15) %>% 
  filter(samp_group < 15) %>%
  filter(treatment != "P-Neg")%>%
  select(treatment, gro_index_t, samp_group,d_s_tran) %>%
  group_by(treatment,samp_group) %>%
  summarise(tot_obs = n(), grow = sum(gro_index_t >= 1)) %>%
  mutate(grow_obs_ind = grow/tot_obs)

ggplot(Tot_grothresh, aes(x=as.factor(samp_group), y = grow_obs_ind)) +
  geom_col()+
  facet_wrap(~treatment, ncol=4) #this plot shows the average viability within each day x treatment group 
