######################
# Project Principles of Statistical data analysis
# 2019 - 2020
# Group 9 : Jan Alexander, Paul Morbee, Joren Verbeke, Steven Wallaert
######################

# Temporary: Load Jan's code to initialize and clean armpit data frame

packages <- c('coin','ggplot2', 'reshape2', 'dplyr', 'gridGraphics')
lapply(packages, library, character.only = TRUE)

setwd("C:/Users/paulm/OneDrive/Mastat/Principles of statistical data analysis/Homeworks/Project 20191220")
filename <- "armpit.txt"
armpit <- read.table(filename)
# Step 1: clean the data. We have 5 genders now
# convert the gender 'F ' to 'F'
armpit[armpit$Gender == 'F ', 'Gender'] = 'F'
armpit[armpit$Gender != 'F', 'Gender'] = 'M'
armpit$Gender <- factor(armpit$Gender)
# drop the observation where gender and age were not filled out
armpit <- armpit[!is.na(armpit$Age),]
armpit <- armpit %>% 
        mutate(Corynebacterium.total = rowSums(.[1:4])) %>%
        mutate(Staphylococcus.total = rowSums(.[5:8])) %>%
        mutate(Bacteria.total = Corynebacterium.total + Staphylococcus.total)
armpit$Agecat <- cut(x = armpit$Age, 
                     breaks = c(0, 25, 35, 50, +Inf), 
                     labels = c("under 25", "25-34", "35 - 49", "50 or older"), 
                     ordered_result = T )
##################
#Characteristics of the categorical data 
##################

ftable(armpit$Gender,armpit$BMI, armpit$Agecat)

#We observe 16 classes where a majority of the classes don't conain sufficient data to draw conclusions
#This means we will have to scope: investigating with Gender and BMI, and dropping them from analysis aftwerwards 
#At least one of the two

ftable(armpit$Gender, armpit$Agecat)
#Even leaving aside BMI we still get classes with insufficient data (older women, younger men)

#What with fewer age categories?
armpit$Agecat <- cut(x = armpit$Age, 
                     breaks = c(0, 40, +Inf), 
                     labels = c("under 40", "40 or older"), 
                     ordered_result = T )
ftable(armpit$Gender,armpit$BMI, armpit$Agecat)
#We still observe classes where a number of the classes don't contain sufficient data

############ 
#Conclusions it will beifficult to draw conclusions on the combination of the three factors
#We restore Agecat category as specifically asked to investigate and use the others 
#if they look promising for certain investigations
###########

armpit$Agecat <- cut(x = armpit$Age, 
                     breaks = c(0, 25, 35, 50, +Inf), 
                     labels = c("under 25", "25-34", "35 - 49", "50 or older"), 
                     ordered_result = T )


#Characteristics of the continuous data
#Histograms
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
hist(armpit$Corynebacterium.1)
hist(armpit$Corynebacterium.2)
hist(armpit$Corynebacterium.3)
hist(armpit$Corynebacterium.4)
#We see in all these histos one dominant component with one or more smaller disconnected portions. 
#Some could maybe be attributed to certain categories, others not. In general not enough data to establish a link 
#with confidence

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
hist(armpit$Staphylococcus.1)
hist(armpit$Staphylococcus.2)
hist(armpit$Staphylococcus.3)
hist(armpit$Staphylococcus.4)
#Same pattern of disconnected zones shows up with the staphilos.

#With these histograms we don't expect much normality around here....
#QQPlots
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
qqplot(x=armpit$Staphylococcus.1,y=rnorm(40,mean=mean(armpit$Staphylococcus.1), sd=sd(armpit$Staphylococcus.1)))
qqplot(x=armpit$Staphylococcus.2,y=rnorm(40,mean=mean(armpit$Staphylococcus.2), sd=sd(armpit$Staphylococcus.2)))
qqplot(x=armpit$Staphylococcus.3,y=rnorm(40,mean=mean(armpit$Staphylococcus.3), sd=sd(armpit$Staphylococcus.3)))
qqplot(x=armpit$Staphylococcus.4,y=rnorm(40,mean=mean(armpit$Staphylococcus.4), sd=sd(armpit$Staphylococcus.4)))

#Even the total
layout(matrix(c(1), 1, 1, byrow = TRUE))
qqplot(x=armpit$Staphylococcus.4,y=rnorm(40,mean=mean(armpit$Staphylococcus.total), sd=sd(armpit$Staphylococcus.total)))

#Same story for the Crynebacteria
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
qqplot(x=armpit$Corynebacterium.1,y=rnorm(40,mean=mean(armpit$Corynebacterium.1), sd=sd(armpit$Corynebacterium.1)))
qqplot(x=armpit$Corynebacterium.2,y=rnorm(40,mean=mean(armpit$Corynebacterium.2), sd=sd(armpit$Corynebacterium.2)))
qqplot(x=armpit$Corynebacterium.3,y=rnorm(40,mean=mean(armpit$Corynebacterium.3), sd=sd(armpit$Corynebacterium.3)))
qqplot(x=armpit$Corynebacterium.4,y=rnorm(40,mean=mean(armpit$Corynebacterium.4), sd=sd(armpit$Corynebacterium.4)))

layout(matrix(c(1), 1, 1, byrow = TRUE))
qqplot(x=armpit$Corynebacterium.4,y=rnorm(40,mean=mean(armpit$Corynebacterium.total), sd=sd(armpit$Corynebacterium.total)))

#None of the distributions can be considered as normal. Long tails.

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
boxplot(armpit$Corynebacterium.1)
boxplot(armpit$Corynebacterium.2)
boxplot(armpit$Corynebacterium.3)
boxplot(armpit$Corynebacterium.4)

#for each bacteria a thin IQR and some "outliers", that when composed give a rather familiar image fot the total
#altough not symmetric
layout(matrix(c(1), 1, 1, byrow = TRUE))
boxplot(armpit$Corynebacterium.total)

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
boxplot(armpit$Staphylococcus.1)
boxplot(armpit$Staphylococcus.2)
boxplot(armpit$Staphylococcus.3)
boxplot(armpit$Staphylococcus.4)#De boxplot voor 
#With the exception of Staphylococcus 1 the same look.
layout(matrix(c(1), 1, 1, byrow = TRUE))
boxplot(armpit$Staphylococcus.total)

#What do we see as boxplots for the different categories?
boxplot(armpit$Staphylococcus.total~armpit$Gender)
#Looking at the median we see there is similarity between the categories

boxplot(armpit$Staphylococcus.total~armpit$BMI)
#Looking at the median we see there is similarity between the categories

boxplot(armpit$Staphylococcus.total~armpit$Agecat)
#Here clearly at higher age tend to become less, altough in the zone "35 - 49" it's lower than "50 or older"
#Due to chance, not having well organized the scale or lack of data ? Essentially each zone has > 5 observations.

#We have to be careful with the interpretation of the boxplots e.g.:
boxplot(armpit$Corynebacterium.1~(armpit$Gender+armpit$BMI))
# e.g. Corynebacterium.1~(armpit$Gender+armpit$BMI) is composed with the following observations for M1:
# 96.17158672, 19.2902319, 15.0572351, 0.142247511, 0, 0, 0


######################
# The Statistics
######################

summary(armpit)


mean_bact <- apply(armpit[,1:8],2,mean)
sd_bact <- apply(armpit[,1:8],2,sd)
mean_bact_M <- apply(armpit[armpit$Gender=="M",1:8],2,mean)
sd_bact_M <- apply(armpit[armpit$Gender=="M",1:8],2,sd)
mean_bact_F <- apply(armpit[armpit$Gender=="F",1:8],2,mean)
sd_bact_F <- apply(armpit[armpit$Gender=="F",1:8],2,sd)
mean_bact_BMI_0 <- apply(armpit[armpit$BMI=="0",1:8],2,mean)
sd_bact_BMI_0 <- apply(armpit[armpit$BMI=="0",1:8],2,sd)
mean_bact_BMI_1 <- apply(armpit[armpit$BMI=="1",1:8],2,mean)
sd_bact_BMI_1 <- apply(armpit[armpit$BMI=="1",1:8],2,sd)


df_mean_sd <- data.frame(mean_bact,sd_bact,mean_bact_M,sd_bact_M,mean_bact_F,sd_bact_F,
                         mean_bact_BMI_0,sd_bact_BMI_0,mean_bact_BMI_1,sd_bact_BMI_1)

#We see that for the different bacteria the sd is 2 to 5 times the mean. Combinations could be 
#investigated to see if for combinations of categoral narrower distributions show up
#For the hypotheses we make confidence levels should be clearly stated

#Let's investigate if the categorical variables have some influence. We opt for Kruskal-Wallis (least conditions) 
#Why: 1. if there is a big change between categories it will be noted
#Why: 2. if there is a correlation between the first and the last categorie it will be noted, if significant
#Gender
kruskal.test(armpit$Corynebacterium.1,armpit$Gender)$p.value
kruskal.test(armpit$Corynebacterium.2,armpit$Gender)$p.value
kruskal.test(armpit$Corynebacterium.3,armpit$Gender)$p.value
kruskal.test(armpit$Corynebacterium.4,armpit$Gender)$p.value
kruskal.test(armpit$Corynebacterium.total,armpit$Gender)$p.value
#Nowhere a significant influence detected

kruskal.test(armpit$Staphylococcus.1,armpit$Gender)$p.value
kruskal.test(armpit$Staphylococcus.2,armpit$Gender)$p.value
kruskal.test(armpit$Staphylococcus.3,armpit$Gender)$p.value
kruskal.test(armpit$Staphylococcus.4,armpit$Gender)$p.value
kruskal.test(armpit$Staphylococcus.total,armpit$Gender)$p.value
#Nowhere a significant influence detected

#BMI 
kruskal.test(armpit$Corynebacterium.1,armpit$BMI)$p.value
kruskal.test(armpit$Corynebacterium.2,armpit$BMI)$p.value
kruskal.test(armpit$Corynebacterium.3,armpit$BMI)$p.value
kruskal.test(armpit$Corynebacterium.4,armpit$BMI)$p.value
kruskal.test(armpit$Corynebacterium.total,armpit$BMI)$p.value
#Nowhere a significant influence detected

kruskal.test(armpit$Staphylococcus.1,armpit$BMI)$p.value
kruskal.test(armpit$Staphylococcus.2,armpit$BMI)$p.value
kruskal.test(armpit$Staphylococcus.3,armpit$BMI)$p.value
kruskal.test(armpit$Staphylococcus.4,armpit$BMI)$p.value
kruskal.test(armpit$Staphylococcus.total,armpit$BMI)$p.value
#Nowhere a significant influence detected

#AgeCat
kruskal.test(armpit$Corynebacterium.1,armpit$Agecat)$p.value
kruskal.test(armpit$Corynebacterium.2,armpit$Agecat)$p.value
kruskal.test(armpit$Corynebacterium.3,armpit$Agecat)$p.value
kruskal.test(armpit$Corynebacterium.4,armpit$Agecat)$p.value
kruskal.test(armpit$Corynebacterium.total,armpit$Agecat)$p.value

#Relative abundance the individuala Corynebacteria are not significantly influenced by Age.
#The total is clearly linked to Agecat
#Question: do they work in the same direction? This we work out further with correlations

kruskal.test(armpit$Staphylococcus.1,armpit$Agecat)$p.value
kruskal.test(armpit$Staphylococcus.2,armpit$Agecat)$p.value
kruskal.test(armpit$Staphylococcus.3,armpit$Agecat)$p.value
kruskal.test(armpit$Staphylococcus.4,armpit$Agecat)$p.value
kruskal.test(armpit$Staphylococcus.total,armpit$Agecat)$p.value

#The same pattern (but with opposite direction pops up for the the Staphylos)


#Do they work in the same direction: To find out we need to use correlation

cor.test(armpit$Corynebacterium.total,armpit$Age,method = "pearson")$estimate
cor.test(armpit$Corynebacterium.total,armpit$Age,method = "kendall")$estimate
cor.test(armpit$Corynebacterium.total,armpit$Age,method = "spearman")$estimate
#All are positively correlated; 0.47, 0.26, 0.39, with pearson having a 95% confidence interval without 0 in it

#What can we say about the individual components?
#1
cor.test(armpit$Corynebacterium.1,armpit$Age,method = "pearson")$estimate
cor.test(armpit$Corynebacterium.1,armpit$Age,method = "kendall")$estimate
cor.test(armpit$Corynebacterium.1,armpit$Age,method = "spearman")$estimate
#All are positively correlated; 0.28, 0.15, 0.21, with pearson having a 95% confidence interval with 0 in it

#2
cor.test(armpit$Corynebacterium.2,armpit$Age,method = "pearson")$estimate
cor.test(armpit$Corynebacterium.2,armpit$Age,method = "kendall")$estimate
cor.test(armpit$Corynebacterium.2,armpit$Age,method = "spearman")$estimate
#All are positively correlated; 0.36, 0.19, 0.28, with pearson having a 95% confidence interval without 0 in it

#3
cor.test(armpit$Corynebacterium.3,armpit$Age,method = "pearson")$estimate
cor.test(armpit$Corynebacterium.3,armpit$Age,method = "kendall")$estimate
cor.test(armpit$Corynebacterium.3,armpit$Age,method = "spearman")$estimate
#Correlations; -0.11, 0.18, 0.23, with pearson having a 95% confidence interval with 0 in it

#4
cor.test(armpit$Corynebacterium.4,armpit$Age,method = "pearson")$estimate
cor.test(armpit$Corynebacterium.4,armpit$Age,method = "kendall")$estimate
cor.test(armpit$Corynebacterium.4,armpit$Age,method = "spearman")$estimate

#Corretalations; 0.36, 0.17, 0.23, with pearson having a 95% confidence interval without 0 in it
#With the excetion Pearson (we will discuss further on) all are relates with age and enforce each other

#What about the AgeCat
#total
cor.test(armpit$Corynebacterium.total,as.numeric(armpit$Agecat),method = "pearson")$estimate
cor.test(armpit$Corynebacterium.total,as.numeric(armpit$Agecat),method = "kendall")$estimate
cor.test(armpit$Corynebacterium.total,as.numeric(armpit$Agecat),method = "spearman")$estimate
#All positively gecorrelated, 0.50, 0.32, 0.43 with a lot of ties for Spearman and Kendall

#1
cor.test(armpit$Corynebacterium.1,as.numeric(armpit$Agecat),method = "pearson")$estimate
cor.test(armpit$Corynebacterium.1,as.numeric(armpit$Agecat),method = "kendall")$estimate
cor.test(armpit$Corynebacterium.1,as.numeric(armpit$Agecat),method = "spearman")$estimate
#All positively (weakly) gecorrelated, 0.27, 0.20, 0.25 and for Pearson 0 in the 95% confidence interval

#2
cor.test(armpit$Corynebacterium.2,as.numeric(armpit$Agecat),method = "pearson")$estimate
cor.test(armpit$Corynebacterium.2,as.numeric(armpit$Agecat),method = "kendall")$estimate
cor.test(armpit$Corynebacterium.2,as.numeric(armpit$Agecat),method = "spearman")$estimate
#All positively (moderately) gecorrelated, 0.41, 0.23, 0.30 and for Pearson 0 not in the 95% confidence interval
#This is the strongest component in the correlation

#3
cor.test(armpit$Corynebacterium.3,as.numeric(armpit$Agecat),method = "pearson")$estimate
cor.test(armpit$Corynebacterium.3,as.numeric(armpit$Agecat),method = "kendall")$estimate
cor.test(armpit$Corynebacterium.3,as.numeric(armpit$Agecat),method = "spearman")$estimate
#Kendall, Spearman (weakly) gecorrelated, -0.10, 0.20, 0.25 and for Pearson 0 in the 95% confidence interval

#4
cor.test(armpit$Corynebacterium.4,as.numeric(armpit$Agecat),method = "pearson")$estimate
cor.test(armpit$Corynebacterium.4,as.numeric(armpit$Agecat),method = "kendall")$estimate
cor.test(armpit$Corynebacterium.4,as.numeric(armpit$Agecat),method = "spearman")$estimate
#Kendall, Spearman (moderately) gecorrelated, 0.31, 0.17, 0.22 and for Pearson 0 not in the 95% confidence interval



#Are we entitled to use Pearson with Agecat ?
scatter.smooth(as.numeric(armpit$Agecat),armpit$Corynebacterium.total)
scatter.smooth(as.numeric(armpit$Agecat),armpit$Corynebacterium.3) #Component had negative Pearson

#Are we entitled to use Pearson with Age ?
scatter.smooth(armpit$Age,armpit$Corynebacterium.total)
scatter.smooth(armpit$Age,armpit$Corynebacterium.3)

#####################################
#What still has to be investigated: confidence intervals for the correlation coefficients.
##95% confidence intervals Spearman and Kendall ? How?
#In progress: I am looking at bootstrapping for the moment
#####################################