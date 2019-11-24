######################
# Project Principles of Statistical data analysis
# 2019 - 2020
# Group 9 : Jan Alexander, Paul Morbée, Joren Verbeke, Steven Wallaert
######################

#Part 1: data cleaning and initial data analysis

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

# Step 2: Make new variables: total of all Corynebacterium species, total of all Staphylococcus species and total of all bacteria
armpit <- armpit %>% 
  mutate(Corynebacterium.total = rowSums(.[1:4])) %>%
  mutate(Staphylococcus.total = rowSums(.[5:8])) %>%
  mutate(Bacteria.total = Corynebacterium.total + Staphylococcus.total)
summary(armpit)

# Step 3: Distributions of the continious variables

AgeDist <- ggplot(data = armpit, aes(x = Age)) + 
  geom_histogram(color = 'black',
                 fill = 'blue') +
  ggtitle("Histogram of patient age") +
  labs(x='Patient age [years]', 
       y='Frequency')
AgeDist

Corynebacterium.totalDist <- ggplot(data = armpit, aes(x = Corynebacterium.total)) + 
  geom_histogram(color = 'black',
                 fill = 'blue',
                 bins = 10) +
  ggtitle("Histogram of relative abundance of Corynebacterium species") +
  labs(x='Corynebacterium spp. [%]', 
       y='Frequency')
Corynebacterium.totalDist

#Corynebacterium relative abundance does not follow a normal distribution.  The genus was absent in most subjects.
#No histogram was made for Staphylococci as this equals 100-Corynebacterium.total


#Part 2 genus composition change with age, strategy 1

#Protocol: a new variable was made with 4 age categories. 
# A boxplot and summary table was made of the relative abundance of Corynebacterium per age category.
# Finally, a Kruskal-Wallis test was performed to study if higher values were more likely in certain groups.

# New variable: Agecat --> 4 age categories defined
armpit$Agecat <- cut(x = armpit$Age, 
                     breaks = c(0, 25, 35, 50, +Inf), 
                     labels = c("under 25", "25-34", "35 - 49", "50 or older"), 
                     ordered_result = T )

# Distributions of the subjects over these age categories
AgeCatDist <- ggplot(armpit, aes(x = Agecat)) +
  geom_bar() + 
  ggtitle('Number of observations in each age category') +
  labs(x='Age category',
       y='Number of observations')
AgeCatDist

# Boxplots of relative abundance of Corynebacterium per age category.
boxcor <- ggplot(armpit, aes(x = Agecat, y=Corynebacterium.total ))+
  geom_boxplot() +
  ggtitle('boxplots of corynebacterium abundance per age category') +
  labs(x='Age category',
       y='Relative abundance of Corynebacterium (%)')

boxcor

#With 2 exceptions, Corynebacterium spp. abundance low in young age categories. In subjects older than 35 years there was a lot of variation. 

# Summary table per age category

tablecor <- armpit %>%
  group_by(Agecat) %>%
  summarize(n = n(),
            mean = mean(Corynebacterium.total, na.rm = TRUE),
            min = min(Corynebacterium.total, na.rm = TRUE),
            median = median(Corynebacterium.total, na.rm = TRUE),
            max = max(Corynebacterium.total, na.rm = TRUE),
            n_positive = sum(Corynebacterium.total>0, na.rm = TRUE),
            percentage_positive=n_positive/n*100)
tablecor

#Kruskal-Wallis test: test if higher values are more likely in certain groups 

kruskal.test(Corynebacterium.total ~ Agecat, data = armpit)

#Part 2 genus composition change with age, strategy 2

#Protocol: a scatterplot with loess curve was made of the relative abundance of Corynebacterium by age. 
#Pearson and Spearman correlation coefficients were calculated and tested.

# Scatterplot to get a first idea how relative abundance of Corynebacterium evolves with age
corbyage <- ggplot(armpit, aes(x=Age, y=Corynebacterium.total)) + 
  geom_point()+
  geom_smooth() +
  ggtitle('relative abundance of Corynebacterium vs subject age') +
  labs(x = 'Age [years]', 
       y = 'Relative abundance of Corynebacterium') 
corbyage


# What does the correlation function in R actually deliver as an output?
cor(armpit$Age, armpit$Corynebacterium.total, method = "pearson")
cor(armpit$Age, armpit$Corynebacterium.total, method = "spearman")

cor.test(armpit$Age, armpit$Corynebacterium.total, 
         alternative = c('two.sided'),
         method = 'pearson',
         conf.level = 0.95)

cor.test(armpit$Age, armpit$Corynebacterium.total, 
         alternative = c('two.sided'),
         method = 'kendall',
         conf.level = 0.95)

#TBA: (dis)advantages of strategy 1 and 2
