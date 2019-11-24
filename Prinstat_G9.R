######################
# Project Principles of Statistical data analysis
# 2019 - 2020
# Group 9 : Jan Alexander, Paul Morbée, Joren Verbeke, Steven Wallaert
######################

# Load relevant packages:
packages <- c('coin','ggplot2', 'reshape2', 'dplyr', 'gridGraphics')
lapply(packages, library, character.only = TRUE)

# Load the file 'armpit.txt', containing the data
filename <- "armpit.txt"
armpit <- read.table(filename)
summary(armpit)

# Step 1: assure we have correct gender labels, now we have 5.
# convert the gender 'F ' to 'F' and ' M' to 'M'
armpit[grepl('F', armpit$Gender), 'Gender'] <- 'F'
armpit[grepl('M', armpit$Gender), 'Gender'] <- 'M'
# drop the observation where gender and age were not filled out --> TODO : is this acceptable?
armpit <- armpit[armpit$Gender != '',]
armpit$Gender <- factor(armpit$Gender)
armpit <- armpit[!is.na(armpit$Age),]
summary(armpit)

# Make new variables: total of all Corynebacterium species
# Make new variables: total of all Staphylococcus species
armpit <- armpit %>% 
  mutate(Corynebacterium.total = rowSums(.[1:4])) %>%
  mutate(Staphylococcus.total = rowSums(.[5:8])) %>%
  mutate(Bacteria.total = Corynebacterium.total + Staphylococcus.total)
summary(armpit)

# Distributions of the patient age
AgeDist <- ggplot(data = armpit, aes(x = Age)) + 
  geom_histogram(color = 'black',
                 fill = 'blue') +
  ggtitle("Histogram of patient age") +
  labs(x='Patient age [years]', 
       y='Frequency')
AgeDist

# New variable: Agecat --> 4 age categories defined
armpit$Agecat <- cut(x = armpit$Age, 
                     breaks = c(0, 25, 35, 50, +Inf), 
                     labels = c("under 25", "25-34", "35 - 49", "50 or older"), 
                     ordered_result = T )

# Distributions of the patients over these age categories
AgeCatDist <- ggplot(armpit, aes(x = Agecat)) +
  geom_bar() + 
  ggtitle('Number of observations in each age category') +
  labs(x='Age category',
       y='Number of observations')
AgeCatDist

# Scatterplot to get a first idea how concentrations of both bacteria genera evolve with age
AgeEvol <- ggplot(armpit, aes(x = Age)) +
  geom_point(aes(y = Staphylococcus.total), color = 'blue') +
  geom_point(aes(y = Corynebacterium.total), color = 'red') + 
  ggtitle('All observations: total of Staphylococcus and 
          total of Corynebacterium vs patient age') +
  labs(x = 'Age [years]', 
       y = 'genus proportion') 
  
# Boxplots with the same objectives as the scatterplot for the age categories.
boxstaph <- ggplot(armpit, aes(x = Agecat, y=Staphylococcus.total ))+
  geom_boxplot() +
  ggtitle('boxplots of staphylococcus proportion vs age category') +
  labs(x='Age category', 
       y = 'Proportion of Staphylococcus')

boxstaphGender <- ggplot(armpit, aes(x = Gender, y=Staphylococcus.total ))+
  geom_boxplot() +
  ggtitle('boxplots of staphylococcus proportion vs gender category') +
  labs(x='Gender', 
       y = 'Proportion of Staphylococcus')

boxcor <- ggplot(armpit, aes(x = Agecat, y=Corynebacterium.total ))+
  geom_boxplot() +
  ggtitle('boxplots of corynebacterium proportion vs age category') +
  labs(x='Age category',
       y='Proportion of Corynebacterium')

AgeEvol
boxstaph
boxstaphGender
boxcor

# What does the correlation function in R actually deliver as an output?
cor(armpit$Age, armpit$Staphylococcus.total, method = "pearson")
cor(armpit$Age, armpit$Staphylococcus.total, method = "spearman")

cor.test(armpit$Age, armpit$Staphylococcus.total, 
         alternative = c('two.sided'),
         method = 'pearson',
         conf.level = 0.95)

cor.test(armpit$Agecat, armpit$Staphylococcus.total, 
         alternative = c('two.sided'),
         method = 'kendall',
         conf.level = 0.95)
