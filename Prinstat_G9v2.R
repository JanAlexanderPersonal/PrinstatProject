######################
# Project Principles of Statistical data analysis
# 2019 - 2020
# Group 9 : Jan Alexander, Paul Morbée, Joren Verbeke, Steven Wallaert
######################

figure.width <- 3.2
figure.height <- 3

#Part 1: data cleaning and initial data analysis

# Load relevant packages:
packages <- c('coin','ggplot2', 'reshape2', 'dplyr', 'gridGraphics', 'tidyverse', 'corrplot', 'tidyr', 'tikzDevice', 'stargazer')
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
armpit <- armpit[!is.na(armpit$Age),]
armpit$Gender <- factor(armpit$Gender)
armpit$BMI <- factor(armpit$BMI)
str(armpit)
# From this summary, we can see that the bacteria total is 100% for all observations
summary(armpit)

# Step 2: Make new variables: total of all Corynebacterium species, total of all Staphylococcus species and total of all bacteria
armpit <- armpit %>% 
  mutate(Corynebacterium.total = rowSums(.[1:4])) %>%
  mutate(Staphylococcus.total = rowSums(.[5:8])) %>%
  mutate(Bacteria.total = Corynebacterium.total + Staphylococcus.total)
summary(armpit)



# Step 2.2: New variable Age above 40 or not

armpit <- armpit %>%
  mutate(age.40plus = factor(ifelse(Age > 40, 1, 0)))

# Step 3: Distributions of the continious variables
AgeDist <- ggplot(data = armpit, aes(x = Age)) + 
  geom_histogram(color = 'black',
                 fill = 'blue',
                 binwidth = 5) +
  ggtitle("Histogram of patient age") +
  labs(x='Patient age [years]', 
       y='Frequency')

tikz(file = 'plot_AgeDistribution.tex', standAlone = FALSE, width = figure.width, height = figure.height)
AgeDist
dev.off()

Corynebacterium.totalDist <- ggplot(data = armpit, aes(x = Corynebacterium.total)) + 
  geom_histogram(color = 'black',
                 fill = 'blue',
                 bins = 10) +
  ggtitle("Histogram of relative abundance of Corynebacterium species") +
  labs(x='Corynebacterium spp. [\\%]', 
       y='Frequency')
tikz(file = 'plot_Coryne_totalDist.tex', standAlone = FALSE, width = figure.width, height = figure.height)
Corynebacterium.totalDist
dev.off()

#Corynebacterium relative abundance does not follow a normal distribution.  The genus was absent in most subjects.
#No histogram was made for Staphylococci as this equals 100-Corynebacterium.total

wilcox.test(Corynebacterium.total ~ BMI, data = armpit)
wilcox.test(Corynebacterium.total ~ Gender, data = armpit)

#Part 2 genus composition change with age, strategy 1

#Protocol: a new variable was made with 4 age categories. 
# A boxplot and summary table was made of the relative abundance of Corynebacterium per age category.
# Finally, a Kruskal-Wallis test was performed to study if higher values were more likely in certain groups.

# New variable: Agecat --> 4 age categories defined
armpit$Agecat <- cut(x = armpit$Age, 
                     breaks = c(0, 25, 35, 50, +Inf), 
                     labels = c("under 25", "25-34", "35 - 49", "50 or older"), 
                     ordered_result = T )

# Tables to investigate the joint dependencies
tab_BMI_Age <- armpit %>%
  group_by(BMI, Agecat) %>%
  summarise(Corynebacterium.AvgAbundance = mean(Corynebacterium.total), 
            Observations = n()) 
tab_Gender_Age <- armpit %>%
  group_by(Gender, Agecat) %>%
  summarise(Corynebacterium.AvgAbundance  = mean(Corynebacterium.total), 
            Observations = n()) 

plot_BMI_Age <- tab_BMI_Age %>%
  ggplot(aes(x = Agecat, y=Corynebacterium.AvgAbundance , size = Observations, color = BMI)) + 
  geom_point(alpha = 0.75) + 
  theme(legend.position = "bottom") +
  scale_size_continuous(limits = c(1, 10), breaks=c(1, 5, 10))+
  ylim(0, 60)+
  ylab('Average relative abundance \n of Corynebacterium genus [\\%]')+
  xlab('Age category')+
  ggtitle('Relative abundance of Corynebacterium \n as a function of age category')


plot_Gender_age <- tab_Gender_Age %>%
  ggplot(aes(x = Agecat, y=Corynebacterium.AvgAbundance , size = Observations, color = Gender)) + 
  geom_point(alpha = 0.75) + 
  theme(legend.position = "bottom") +
  ylim(0, 60)+
  scale_size_continuous(limits = c(1, 10), breaks=c(1, 5, 10))+
  ylab('Average relative abundance \n of Corynebacterium genus [\\%]')+
  xlab('Age category')+
  ggtitle('Relative abundance of Corynebacterium \n as a function of age category')

tikz(file = 'plot_BMI_Age.tex', standAlone = FALSE, width = figure.width, height = figure.height)
plot_BMI_Age
dev.off()
tikz(file = 'plot_Gender_age.tex', standAlone = FALSE, width = figure.width, height = figure.height)
plot_Gender_age
dev.off()

# Distributions of the subjects over these age categories
AgeCatDist <- ggplot(armpit, aes(x = Agecat)) +
  geom_bar() + 
  ggtitle('Number of observations in each age category') +
  labs(x='Age category',
       y='Number of observations')

tikz(file = 'plot_AgeCat_distribution.tex', standAlone = FALSE, width = figure.width, height = figure.height)
AgeCatDist
dev.off()

# Boxplots of relative abundance of Corynebacterium per age category.
boxcorr_cory_agecat <- ggplot(armpit, aes(x = Agecat, y=Corynebacterium.total ))+
  geom_boxplot() +
  ggtitle('boxplots of corynebacterium abundance per age category') +
  labs(x='Age category',
       y='Relative abundance of Corynebacterium (\\%)')

tikz(file = 'plot_boxplotcorr_Coryne_AgeCat', standAlone = FALSE, width = figure.width, height = figure.height)
boxcorr_cory_agecat
dev.off()

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
stargazer(tablecor, summary = FALSE, out='table_cor.tex', align = TRUE, digits = 2)

#Kruskal-Wallis test: test if higher values are more likely in certain groups 

kruskal.test(Corynebacterium.total ~ Agecat, data = armpit)

#Part 2 genus composition change with age, strategy 2

#Protocol: a scatterplot with loess curve was made of the relative abundance of Corynebacterium by age. 
#Pearson and Spearman correlation coefficients were calculated and tested.

# Scatterplot to get a first idea how relative abundance of Corynebacterium evolves with age
corbyage <- ggplot(armpit, aes(x=Age, y=Corynebacterium.total)) + 
  geom_point()+
  geom_smooth() + # Ziet er mooi uit, maar is dit ook 'correct' (aanvaardbaar)
  ggtitle('relative abundance of Corynebacterium vs subject age') +
  labs(x = 'Age [years]', 
       y = 'Relative abundance of Corynebacterium') 
tikz(file = 'plot_scatterplot_Age_corby.tex', standAlone = FALSE, width = figure.width, height = figure.height)
corbyage
dev.off()


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

# Part 3 relative abundances of the 8 species

#First step: look how many species occur together

armpit$nspecies <-  rowSums(armpit[,c(1:8)]>0)
summary(armpit)

nspeciesDist <- ggplot(data = armpit, aes(x = nspecies)) + 
  geom_histogram(color = 'black',
                 fill = 'blue',
                 bins=7) +
  ggtitle("Histogram of number of concurrent species") +
  labs(x='number of species', 
       y='Frequency')

tikz(file = 'plot_JointOccurance.tex', standAlone = FALSE, width = figure.width, height = figure.height)
nspeciesDist
dev.off()

#In one subject, only 1 species was found.

#Prevalence of each species is shown in the table below

armpit$subjectID <- rownames(armpit)
long_armpit <- armpit %>% gather(Species, Abundance, Corynebacterium.1:Staphylococcus.4)
head(long_armpit,50)

tableprevalence <- long_armpit %>%
  group_by(Species) %>%
  summarize(n_positive = sum(Abundance>0, na.rm = TRUE),
            percentage_positive=n_positive/length(Abundance)*100,
            min_abundance = min(Abundance, na.rm = TRUE),
            max_abundance = max(Abundance, na.rm = TRUE))

tableprevalence


#The table below shows the conditional probabilities of detecting a species when another species was detected

mean(armpit$Corynebacterium.1[armpit$Corynebacterium.2>0]>0)

#TBA: how can we do this for each species vs each other species?

#The plot below shows the spearman correlation coefficients

speciescor <- cor(armpit[ ,c(1:8)],method="spearman")
round(speciescor,2)
tikz(file = 'plot_correlation', standAlone = FALSE, width = 8, height = 8)
corrplot(speciescor, type = "upper",  
         tl.col = "black", tl.srt = 45)
dev.off()


### steven

# dichotomiseren van Coryne abundance
armpit <- armpit %>%
  mutate(Coryne.dichotoom = factor(ifelse(Corynebacterium.total>50, 1, 0)))

plot_dichotoom <- ggplot(armpit, aes(x=Coryne.dichotoom))+ 
  geom_histogram(color = 'black',
                 fill = 'blue') +
  ggtitle("dichotoom") 
tikz(file = 'plot_dichotoom.tex', height = figure.height, width = figure.width, standAlone = FALSE)
plot_dichotoom
dev.off()


theme_set(theme_light())

plot_1 <- ggplot(armpit,aes(x = Age, y = Corynebacterium.total)) +
  geom_point() +
  geom_smooth()

tikz(file = 'plot1.tex', height = figure.height, width = figure.width, standAlone = FALSE)
plot_1
dev.off()

cor.test(armpit$Age, armpit$Corynebacterium.total,
         alternative = "two.sided",
         method = "kendall",
         conf.level = .95)

plot.2 <- ggplot(armpit,aes(x = Age, y = Coryne.dichotoom)) +
  geom_point() +
  geom_smooth()

tikz(file = 'plot2.tex', height = figure.height, width = figure.width, standAlone = FALSE)
plot_2
dev.off()

# dit is niet OK
# cor.test(armpit$Age, armpit$Coryne.dichotoom,
#         alternative = "two.sided",
#         method = "pearson",
#         conf.level = .95)



# install.packages("polycor")
# library(polycor)
# Hm, warning: based on the assumption that the joint distribution 
# of the quantitative variable and a latent continuous variable underlying
# the ordinal variable is bivariate normal
# polyserial(armpit$Age, armpit$Coryne.dichotoom)


### Discrete

cor.test(as.numeric(armpit$Agecat), armpit$Corynebacterium.total,
         alternative = "two.sided",
         method = "kendall",
         conf.level = .95)

tikz(file = 'plot3.tex', height = figure.height, width = figure.width, standAlone = FALSE)
armpit %>%
  select(Corynebacterium.total, Staphylococcus.total, Agecat) %>%
  pivot_longer(1:2, names_to = "Bacterium", values_to = "relabundance") %>%
  ggplot(aes(x = Agecat, y=relabundance, fill=Bacterium)) +
  geom_col(position = "fill")
dev.off()

