######################
# Project Principles of Statistical data analysis
# 2019 - 2020
# Group 9 : Jan Alexander, Paul Morbée, Joren Verbeke, Steven Wallaert
######################

figure.width <- 3.2
figure.height <- 3


round_df <- function(x, digits = 2) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

#Part 1: data cleaning and initial data analysis

# Load relevant packages:
packages <- c('coin','ggplot2', 'reshape2', 'dplyr', 'gridGraphics', 'tidyverse', 'corrplot', 'tidyr', 'tikzDevice', 'stargazer')
lapply(packages, library, character.only = TRUE)

theme_set(theme_light())

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
armpit$BMI <- factor(ifelse(armpit$BMI == 0, 'BMI <= 25', 'BMI > 25'), ordered = TRUE)
str(armpit)
# From this summary, we can see that the bacteria total is 100% for all observations
summary(armpit)

# Step 2: Make new variables: total of all Corynebacterium species, total of all Staphylococcus species and total of all bacteria + New variable Age above 40 or not
armpit <- armpit %>% 
  mutate(Corynebacterium.total = rowSums(.[1:4])) %>%
  mutate(Staphylococcus.total = rowSums(.[5:8])) %>%
  mutate(Bacteria.total = Corynebacterium.total + Staphylococcus.total) %>%
  mutate(Agecat = factor(ifelse(Age > 40, 'over 40', '40 or younger')))
summary(armpit)

# Count the number of observations in both gender en BMI categories
tab_BMI_Gender <- armpit %>%
  group_by(BMI, Gender) %>%
  summarise(Corynebacterium.AvgAbundance = mean(Corynebacterium.total), 
            Observations = n()) 

stargazer(tab_BMI_Gender %>% 
            select(BMI, Gender, Observations) %>%
            spread(data = ., key = BMI, value = Observations), 
          summary = FALSE, nobs = TRUE, out='table_Observations_BMI_Gender.tex')

# Explore the independent variable age

AgeDist <- ggplot(data = armpit, aes(x = Age)) + 
  geom_histogram(color = 'black',
                 fill = 'blue',
                 binwidth = 5) +
  labs(x='Patient age [years]', 
       y='Observation frequency')

stargazer(armpit %>% select(Age) %>% summarize(min = min(Age, na.rm = TRUE),
                                               Q1 = quantile(Age, .25),
                                               median = median(Age, na.rm = TRUE),
                                               mean = mean(Age, na.rm = TRUE),
                                               Q3 = quantile(Age, .75),
                                               max = max(Age, na.rm = TRUE))
          , summary = FALSE, out = 'table_Age_statistics.tex', rownames = FALSE)
tikz(file = 'plot_AgeDistribution.tex', standAlone = FALSE, width = figure.width, height = figure.height)
AgeDist
dev.off()

#Mainly young subjects participated in the study. The age distribution was not normal.

# Step 3: Explore the variables of interest

#Basic measures of location and spread are shown below.

armpit$subjectID <- rownames(armpit)
species_armpit <- armpit %>% gather(Species, Abundance, Corynebacterium.1:Staphylococcus.4) %>%
  mutate(Genus = ifelse(grepl("Cory", Species), "Corynebacterium", "Staphyloccus"))
head(species_armpit,50)


measures_species <- species_armpit %>%
  group_by(Species) %>%
  summarize(mean = mean(Abundance, na.rm = TRUE),
            median = median(Abundance, na.rm = TRUE),
            sd = sd(Abundance, na.rm = TRUE),
            min = min(Abundance, na.rm = TRUE),
            max = max(Abundance, na.rm = TRUE))

genus_armpit <- armpit %>% gather(Species, Abundance, Corynebacterium.total:Staphylococcus.total)
head(genus_armpit,50)

measures_genus <- genus_armpit %>%
  group_by(Species) %>%
  summarize(mean = mean(Abundance, na.rm = TRUE),
            median = median(Abundance, na.rm = TRUE),
            sd = sd(Abundance, na.rm = TRUE),
            min = min(Abundance, na.rm = TRUE),
            max = max(Abundance, na.rm = TRUE))
stargazer(round_df(rbind(measures_genus,measures_species), 1), summary = FALSE, nobs = FALSE, out='table_Statistics_Bacteria.tex', align = FALSE, rownames = FALSE)


#The large differences between means and medians indicate that the values are not normally distributed. This was further studied by plotting histograms.

Corynebacterium.totalDist <- ggplot(data = armpit, aes(x = Corynebacterium.total)) + 
  geom_histogram(color = 'black',
                 fill = 'blue',
                 bins = 10) +
  labs(x='Corynebacterium genus relative abundance [\\%]', 
       y='Frequency')

tikz(file = 'plot_Coryne_totalDist.tex', standAlone = FALSE, width = figure.width, height = figure.height)
Corynebacterium.totalDist
dev.off()

# TODO: What do we prefer: boxplot or histogram?

Corynebacterium.totalDist <- ggplot(data = armpit, aes(y = Corynebacterium.total)) + 
  geom_boxplot(notch = FALSE) + xlab('Relative abundance [\\%]')+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  coord_flip()
  
tikz('boxplot_Cory.tex', width = 4, height = 1.2)
Corynebacterium.totalDist
dev.off()

#Corynebacterium relative abundance does not follow a normal distribution.  The genus was absent in most subjects.
#No histogram was made for Staphylococci as this equals 100-Corynebacterium.total

#Boxplots and individual observations of the relative abundance were plotted per species.

SpeciesDist <- ggplot(species_armpit, aes(Species, Abundance)) +
  geom_boxplot(colour = "grey50",outlier.shape=NA) +
  geom_jitter()+
  geom_line(data=measures_species, aes(x=Species, y=mean,group = 1),color='blue') + 
  geom_line(data=measures_species, aes(x=Species, y=median,group = 1),color='red') + 
  ylab("Relative abundance (%)")+
  xlab("Species")+
  theme(axis.text.x = element_text(angle = 45))+
  annotate(geom="text", x=7, y=80, label="median", color="red")+
  annotate(geom="text", x=7, y=60, label="mean", color="blue")

SpeciesDist



ggplot(species_armpit, aes(Abundance)) +
  geom_histogram(bins = 20) +
  geom_vline(data = measures_species, aes(xintercept = mean), color="blue") +
  geom_vline(data = measures_species, aes(xintercept = median), color = "red") +
  facet_wrap(~Species, nrow = 2, scales = "free") +
  labs(x="Relative species abundance (%)",
       y="Count")

ggplot(species_armpit, aes(x = Species, y = Abundance)) +
  geom_boxplot(aes(color = Genus)) +
  geom_point(data = measures_species, aes(x = Species, y=`mean`), shape = 6, size = 3) +
  coord_flip() +
  labs(x="Species",
       y="Relative species abundance (%)") +
  guides(color = guide_legend(reverse = TRUE), shape = guide_legend(title = "Mean")) 

  
#The figure shows that Staphylococcus 1 was the most common species. Other species were often absent in subjects.
#For each species, abundance did not follow a normal distribution. Mean and median appear to be poor measures of location. 

#Part 2 genus composition change with age, strategy 1

#Protocol: a new variable was made with 4 age categories. 
# A boxplot and summary table was made of the relative abundance of Corynebacterium per age category.
# Finally, a Kruskal-Wallis test was performed to study if higher values were more likely in certain groups.

# Tables to investigate the joint dependencies --> better with combined box plots?
stargazer(armpit %>%
  group_by(BMI, Gender,Agecat) %>%
  summarise(Observations = n())  %>% 
  spread(data = ., key = Agecat, value = Observations),
summary = FALSE, nobs = TRUE, rownames = FALSE, out='table_DistributionAfterAgeDiscr.tex')

#With 4 exceptions, Corynebacterium spp. abundance low in young age categories. In subjects older than 40 years there was a lot of variation. 
plot_Age <-armpit %>% ggplot(aes(x=Agecat, y=Corynebacterium.total)) + 
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(), aes(shape = Agecat)) + 
  theme(legend.position = "bottom") + 
  ylab('Relative abundance \n of Corynebacterium genus [\\%]') + 
  xlab('Age category')

# Summary table per age category
tableAgeCat <- armpit %>%
  group_by(Agecat) %>%
  summarize(n = n(),
            mean = round(mean(Corynebacterium.total, na.rm = TRUE),1),
            min = round(min(Corynebacterium.total, na.rm = TRUE),1),
            median = round(median(Corynebacterium.total, na.rm = TRUE),1),
            max = round(max(Corynebacterium.total, na.rm = TRUE),1),
            n_present = sum(Corynebacterium.total>0, na.rm = TRUE),
            percentage_present=round(n_present/n*100, 1)) 
stargazer(tableAgeCat, summary = FALSE, out='table_AgeCat_Cory.tex', rownames = FALSE, nobs = FALSE)

#Wilcoxon test: test if the distribution of Corynebacterium abundance differed between young and old subjects. 

wilcox.test(Corynebacterium.total ~ Agecat, data = armpit, alternative = "two.sided")
# The test was highly significant. Higher Corynebacterium abundance was more likely in older subjects.

# How comparable are the two age categories apart from the age difference? 

plot_BMI_Age <-armpit %>% ggplot(aes(x=Agecat, y=Corynebacterium.total, fill = BMI)) + 
  geom_boxplot(outlier.alpha = 0)+
  geom_point(aes(shape = BMI), position = position_jitterdodge()) + 
  theme(legend.position = "bottom") + 
  ylab('Relative abundance \n of Corynebacterium genus [\\%]') + 
  xlab('Age category') + scale_fill_manual(values=c("white", "gray"))

plot_Gender_age <-armpit %>% ggplot(aes(x=Agecat, y=Corynebacterium.total, fill = Gender)) + 
  geom_boxplot(outlier.alpha = 0)+
  geom_point(aes(shape = Gender), position = position_jitterdodge()) + 
  theme(legend.position = "bottom") + 
  ylab('Relative abundance \n of Corynebacterium genus [\\%]') + 
  xlab('Age category') + scale_fill_manual(values=c("white", "gray"))

tikz(file = 'plot_Age.tex', standAlone = FALSE, width = figure.width, height = figure.height)
plot_Age
dev.off()

tikz(file = 'plot_BMI_Age.tex', standAlone = FALSE, width = figure.width, height = figure.height)
plot_BMI_Age
dev.off()

tikz(file = 'plot_Gender_age.tex', standAlone = FALSE, width = figure.width, height = figure.height)
plot_Gender_age
dev.off()

#A Kruskal wallis test was performed to test if distribution of Corynebacterium abundance differed between gender and BMI groups.  

armpit %>%
  select(Corynebacterium.total, BMI, Gender) %>%
  mutate(BMI.Gender = factor(ifelse(BMI == 'BMI <= 25', ifelse(Gender == "F", 0, 1), ifelse(Gender == "F", 2, 3)))) %>%
  kruskal.test(Corynebacterium.total ~ BMI.Gender, data = .)

#This test was non-significant. Confounding by BMI and age is unlikely.  


#Part 2 genus composition change with age, strategy 2

#Protocol: a scatterplot with loess curve was made of the relative abundance of Corynebacterium by age. 
#Pearson and Spearman correlation coefficients were calculated and tested.

# # Scatterplot to get a first idea how relative abundance of Corynebacterium evolves with age
# corbyage <- ggplot(armpit, aes(x=Age, y=Corynebacterium.total)) + 
#   geom_point()+
#   geom_smooth() + # Ziet er mooi uit, maar is dit ook 'correct' (aanvaardbaar)
#   ggtitle('relative abundance of Corynebacterium vs subject age') +
#   labs(x = 'Age [years]', 
#        y = 'Relative abundance of Corynebacterium') 
# tikz(file = 'plot_scatterplot_Age_corby.tex', standAlone = FALSE, width = figure.width, height = figure.height)
# corbyage
# dev.off()
# 
# 
# # What does the correlation function in R actually deliver as an output?
# cor(armpit$Age, armpit$Corynebacterium.total, method = "pearson")
# cor(armpit$Age, armpit$Corynebacterium.total, method = "spearman")
# 
# cor.test(armpit$Age, armpit$Corynebacterium.total, 
#          alternative = c('two.sided'),
#          method = 'pearson',
#          conf.level = 0.95)

cor.test(armpit$Age, armpit$Corynebacterium.total, 
         alternative = c('two.sided'),
         method = 'kendall',
         conf.level = 0.95)

source("tau_cor_conf.R")

tau_cor_conf(armpit, "Age", "Corynebacterium.total")

ggplot(armpit, aes(x = Age, y = Corynebacterium.total)) +
  geom_hex()

ggplot(armpit, aes(x = Age, y = Corynebacterium.total)) +
  geom_point(aes(color = Gender, shape= BMI)) 

ggplot(armpit, aes(x = Age, y = Corynebacterium.total)) +
  stat_density_2d(aes(fill=stat(density)), geom = "raster", contour = FALSE) +
  scale_fill_continuous(type = "viridis") +
  theme_minimal() 


#TBA: (dis)advantages of strategy 1 and 2

# Part 3 relative abundances of the 8 species

#The initial data analysis showed that certain species are often absent in subjects.
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

#In one subject, only 1 species was found. None of the subjects tested positive on all 8 species.


#The table below shows the conditional probabilities of detecting a species when another species was detected

mean(armpit$Corynebacterium.1[armpit$Corynebacterium.2>0]>0)

condprob <- matrix(nrow = 8, ncol = 8)

armpit %>%
  select(Corynebacterium.1:Corynebacterium.4, Staphylococcus.1:Staphylococcus.4) -> bacteria



for(i in 1:8){
  for(j in 1:8){
    condprob[i,j] <- mean(bacteria[i][bacteria[j]>0]>0)
  }
}

rownames(condprob) <- colnames(bacteria)
colnames(condprob) <- colnames(bacteria)

round(condprob,  2)

#Can we do the same thing with a chi square test?

chisq.test(armpit$Corynebacterium.1>0,armpit$Corynebacterium.2>0)$p.value

#For loops seems not to work.



Pchisq <- matrix(nrow = 8, ncol = 8)


for(i in 1:8){
  for(j in 1:8){
    Pchisq[i,j] <- chisq.test(bacteria[i]>0,bacteria[j]>0)$p.value
  }
}

rownames(Pchisq) <- colnames(bacteria)
colnames(Pchisq) <- colnames(bacteria)

round(Pchisq, 2)


#The plot below shows the Kendall correlation coefficients

speciescor <- cor(armpit[ ,c(1:8)],method="kendall")
round(speciescor,2)
tikz(file = 'plot_correlation.tex', standAlone = FALSE, width = 8, height = 8)
corrplot(speciescor, type = "upper",  
         tl.col = "black", tl.srt = 45,method="number")
dev.off()




# dichotomiseren van Coryne abundance
armpit <- armpit %>%
  mutate(Coryne.dichotoom = factor(ifelse(Corynebacterium.total>50, 1, 0)))

plot_dichotoom <- ggplot(armpit, aes(x=Coryne.dichotoom))+ 
  geom_bar(color = 'black',
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

plot_2 <- ggplot(armpit,aes(x = Age, y = as.numeric(Coryne.dichotoom)-1)) +
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

