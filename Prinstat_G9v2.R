######################
# Project Principles of Statistical data analysis
# 2019 - 2020
# Group 9 : Jan Alexander, Paul Morbée, Joren Verbeke, Steven Wallaert
######################

figure.width <- 3.2
figure.height <- 2.8

#####
round_df <- function(x, digits = 2) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

source("tau_cor_conf.R")

# Load packages ####

# Load relevant packages:
packages <- c('reshape2', 'gridGraphics', 'tidyverse', 'corrplot', 'tikzDevice', 'stargazer', 'hexbin')
lapply(packages, library, character.only = TRUE)

theme_set(theme_light())

# Load data, data cleaning and descriptive data analysis ####
# Load the file 'armpit.txt', containing the data
filename <- "armpit.txt"
armpit <- read.table(filename)
str(armpit)

# Step 1: assure we have correct gender labels, now we have 5.
# convert the gender 'F ' to 'F' and ' M' to 'M'
armpit[grepl('F', armpit$Gender), 'Gender'] <- 'F'
armpit[grepl('M', armpit$Gender), 'Gender'] <- 'M'
# drop the observation where gender and age were not filled out. 
# Factor the BMI values with their description.
armpit <- armpit[!is.na(armpit$Age),]
armpit$Gender <- factor(armpit$Gender)
armpit$BMI <- factor(ifelse(armpit$BMI == 0, 'BMI <= 25', 'BMI > 25'), ordered = TRUE)
armpit$subjectID <- as.numeric(rownames(armpit))
str(armpit)
# From this summary, we can see that the bacteria total is 100% for all observations
summary(armpit)

# Make new variables: total of all Corynebacterium species, ####
# total of all Staphylococcus species and total of all bacteria
# New variable: age categories --> Age above 40 or not
armpit <- armpit %>% 
  mutate(Corynebacterium.total = rowSums(.[1:4])) %>%
  mutate(Staphylococcus.total = rowSums(.[5:8])) %>%
  mutate(Bacteria.total = Corynebacterium.total + Staphylococcus.total) %>%
  mutate(Agecat = factor(ifelse(Age > 40, 'over 40', '40 or younger')))
summary(armpit)

# Count the number of observations in both gender en BMI categories
# AvgAbundance is added to provide a first view
tab_BMI_Gender <- armpit %>%
  group_by(BMI, Gender) %>%
  summarise(Corynebacterium.AvgAbundance = mean(Corynebacterium.total), 
            Observations = n()) 

# Export table to include into report (only number of observations)
stargazer(tab_BMI_Gender %>% 
            select(BMI, Gender, Observations) %>%
            spread(data = ., key = BMI, value = Observations), 
          summary = FALSE, nobs = TRUE, out='table_Observations_BMI_Gender.tex')

# AgeDist: histogram of age distribution in dataset
AgeDist <- ggplot(data = armpit, aes(x = Age)) + 
  geom_histogram(color = 'black',
                 fill = 'white',
                 binwidth = 5) +
  scale_y_continuous(breaks = seq(0, 10, 2))+
  labs(x='Subject age [years]', 
       y='\\# Subjects') +
  theme(axis.title=element_text(size=14)) 

# Export table with basic statistics age distribution for report
stargazer(armpit %>% select(Age) %>% summarize(min = min(Age, na.rm = TRUE),
                                               Q1 = quantile(Age, .25),
                                               median = median(Age, na.rm = TRUE),
                                               mean = mean(Age, na.rm = TRUE),
                                               Q3 = quantile(Age, .75),
                                               max = max(Age, na.rm = TRUE))
          , summary = FALSE, out = 'table_Age_statistics.tex', rownames = FALSE)
# Export histogram figure Age distribution
tikz(file = 'plot_AgeDistribution.tex', standAlone = FALSE, width = figure.width, height = figure.height*0.70)
AgeDist
dev.off()


#Mainly young subjects participated in the study. 
# The age distribution does not correspond to the distribution of the Belgian population as a whole.

# Explore the variables of interest (Relative abundance of bacteria species) ####
# Basic measures of location and spread are shown below and exported to a table for the report

# Make a long table with the abundance for each of the 8 species
species_armpit <- armpit %>% 
  gather(Species, Abundance, Corynebacterium.1:Staphylococcus.4) %>%
  mutate(Genus = ifelse(grepl("Cory", Species), "Corynebacterium", "Staphyloccus"))
head(species_armpit,3)

# Make table with basic statistics of these relative species abundances
measures_species <- species_armpit %>%
  group_by(Species) %>%
  summarize(mean = mean(Abundance, na.rm = TRUE),
            median = median(Abundance, na.rm = TRUE),
            sd = sd(Abundance, na.rm = TRUE),
            iqr = IQR(Abundance, na.rm = TRUE),
            min = min(Abundance, na.rm = TRUE),
            max = max(Abundance, na.rm = TRUE))

# Make a long table with the abuncance for both genera (like was done for the species)
genus_armpit <- armpit %>% 
  gather(Species, Abundance, Corynebacterium.total:Staphylococcus.total)
head(genus_armpit,3)

# Make a table with the same basic statistics of the genera abundances as was done for the relative species abundances
measures_genus <- genus_armpit %>%
  group_by(Species) %>%
  summarize(mean = mean(Abundance, na.rm = TRUE),
            median = median(Abundance, na.rm = TRUE),
            sd = sd(Abundance, na.rm = TRUE),
            iqr = IQR(Abundance, na.rm = TRUE),
            min = min(Abundance, na.rm = TRUE),
            max = max(Abundance, na.rm = TRUE))
round_df(rbind(measures_genus,measures_species), 1)

# Export both the table with the species basic statistics and the genera basic statistics
stargazer(round_df(rbind(measures_genus,measures_species), 1), summary = FALSE, nobs = FALSE, out='table_Statistics_Bacteria.tex', align = FALSE, rownames = FALSE)

bacteria_armpit <- armpit %>% 
  gather(Bacteria, Abundance, c(Corynebacterium.1:Corynebacterium.4,
                                                            Corynebacterium.total,
                                                            Staphylococcus.1:Staphylococcus.4,
                                                            Staphylococcus.total)) %>%
  mutate(Genus = ifelse(grepl("Cory", Bacteria), "Corynebacterium", "Staphyloccus"))

measures_bacteria <- bacteria_armpit %>%
  group_by(Bacteria) %>%
  summarize(mean = mean(Abundance, na.rm = TRUE),
            median = median(Abundance, na.rm = TRUE),
            sd = sd(Abundance, na.rm = TRUE),
            iqr = IQR(Abundance, na.rm = TRUE),
            min = min(Abundance, na.rm = TRUE),
            max = max(Abundance, na.rm = TRUE))

#The large differences between means and medians indicate that the values are not normally distributed. 
# This was further studied by plotting a boxplot.

# Corynebacterium relative abundance does not follow a normal distribution.  The genus was absent in most subjects.
# No histogram was made for Staphylococci as this equals 100 % -Corynebacterium.total

# Boxplot indicating the distribution of the relative species abundance
  
Boxplot_bacteria <- ggplot(bacteria_armpit, aes(x = Bacteria, y = Abundance)) +
  geom_boxplot(aes(fill = Genus)) +
  scale_fill_manual(values=c("gray", "white")) +
  geom_point(data = measures_bacteria, 
             aes(x = Bacteria, y=`mean`), shape = 6, size = 3) +
  coord_flip() +
  theme(legend.position = "bottom") +
  labs(x="Bacteria",
       y="Relative bacteria abundance (\\%)") +
  guides(color = guide_legend(reverse = FALSE), 
         shape = guide_legend(title = "Mean")) +
  theme(axis.title=element_text(size=16), 
        axis.text.y = element_text(size = 14)) 

tikz(file = 'plot_Boxplot_species.tex', standAlone = FALSE, width = figure.width * 2, height = figure.height)
  Boxplot_bacteria
dev.off()


# The figure shows that Staphylococcus 1 was the most common species. Other species were often absent in subjects (RSA = 0.0 %).
# For each species, abundance did not follow a normal distribution. 
# Mean and median appear to be poor measures of location for these distributions. 

# Relative abundance genus change with age, strategy 1 : Age as a discrete variable ####

#Protocol: a new variable was made with 2 age categories. 
# A boxplot and summary table was made of the relative abundance of Corynebacterium per age category.
# Finally, a Kruskal-Wallis test was performed to study if higher values were more likely in certain groups.

# Tables to investigate the joint dependencies --> better with combined box plots?
stargazer(armpit %>%
  group_by(BMI, Gender,Agecat) %>%
  summarise(Observations = n())  %>% 
  spread(data = ., key = Agecat, value = Observations),
summary = FALSE, nobs = TRUE, rownames = FALSE, out='table_DistributionAfterAgeDiscr.tex')

# With 4 exceptions, Corynebacterium spp. abundance low in young age categories. 
# In subjects older than 40 years there was a lot of variation. 
plot_Age <-armpit %>% ggplot(aes(x=Agecat, y=Corynebacterium.total)) + 
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width = 0.25), aes(shape = Agecat)) + 
  stat_summary(fun.y = mean, aes(shape = Agecat) ,geom="point",color = 'darkgray', size=3, position = position_dodge(width = 0.1)) + 
  theme(legend.position = "none", axis.title=element_text(size=14)) + 
  ylab('Relative abundance \n of Corynebacterium genus [\\%]') + 
  xlab('Age category') + 
  theme(axis.title=element_text(size=16))

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

# It is important to avoid confounding: 
# How comparable are the two age categories apart from the age difference? 
# It could for example be that the investigated sample only contains young men and women over 40 (not the case, just an example).
# Remark : outlier.alpha = 0 makes the outliers in the boxplot invisible. This is only acceptable in this case
# Because these are indicated by the overlayed point plots!
plot_BMI_Age <-armpit %>% ggplot(aes(x=Agecat, y=Corynebacterium.total, fill = BMI)) + 
  geom_boxplot(outlier.alpha = 0)+
  stat_summary(fun.y = mean, geom="point",aes( shape = BMI),color = 'darkgray', size=3, position = position_jitterdodge())+
  geom_point(aes(shape = BMI), position = position_jitterdodge()) + 
  theme(legend.position = "bottom") + 
  ylab('Relative abundance \n of Corynebacterium genus [\\%]') + 
  xlab('Age category') + scale_fill_manual(values=c("white", "lightgray"))

plot_Gender_age <-armpit %>% ggplot(aes(x=Agecat, y=Corynebacterium.total, fill = Gender)) + 
  geom_boxplot(outlier.alpha = 0)+
  stat_summary(fun.y = mean, geom="point",aes( shape = Gender),color = 'darkgray', size=3, position = position_jitterdodge())+
  geom_point(aes(shape = Gender), position = position_jitterdodge()) + 
  theme(legend.position = "bottom") + 
  ylab('Relative abundance \n of Corynebacterium genus [\\%]') + 
  xlab('Age category') + scale_fill_manual(values=c("white", "lightgray"))

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
wilcox.test(Corynebacterium.total ~ BMI, data = armpit, alternative = "two.sided")
wilcox.test(Corynebacterium.total ~ Gender, data = armpit, alternative = "two.sided")

#These tests were non-significant. Confounding by BMI and age is unlikely.  


# Relative abundance genus change with age, strategy 2: Age as a continuous variable ####

# Function to determine the p-value of the correlation with a permutation test
Kendall.perm.test <- function (x, y){
  nsim <- 10000
  res <- numeric(nsim) #Empty vector with length of simulations
  for (i in 1:nsim) {
    perm <- sample(x)  #Permutation
    res[i] <- cor(perm, y, method = 'kendall') #Correlation in permutated dataset
  }
  obs <-  cor(x, y, method = 'kendall') #Observed correlation in original dataset
  mean(abs(res)>=abs(obs)) #Proportion of permutations with equal or more extreme difference than the observed difference 
}

# Association between age and relative abundance of the Corynebacterium is tested with the Kendall correlation coefficient.

cor.test(armpit$Age, armpit$Corynebacterium.total, 
         alternative = c('two.sided'),
         method = 'kendall',
         conf.level = 0.95)

cor(armpit$Age, armpit$Corynebacterium.total, method = 'kendall')

Kendall.perm.test(armpit$Age, armpit$Corynebacterium.total)

tau_cor_conf(armpit, "Age", "Corynebacterium.total")

plot_scatterplot_age_cory <- ggplot(armpit, aes(x = Age, y = Corynebacterium.total)) +
  geom_point(aes(color = Gender, shape= BMI), size = 3) +
  theme(legend.position = "bottom") + 
  ylab('Relative abundance \n of Corynebacterium genus [\\%]') + 
  xlab('Age [years]') + 
  scale_color_manual(values=c("darkgray", "black"))

tikz(file = 'plot_scatterplot_age_cory.tex', standAlone = FALSE, width = figure.width*2, height = figure.height*2)
plot_scatterplot_age_cory
dev.off()

# Part 3 relative abundances of the 8 species ####

#The initial data analysis showed that certain species are often absent in subjects.
#First step: look how many species occur together

armpit$nspecies <-  rowSums(armpit[,c(1:8)]>0)

nspeciesDist <- ggplot(data = armpit, aes(x = nspecies)) + 
  geom_histogram(color = 'black',
                 fill = 'white',
                 bins=7) +
  scale_x_continuous(breaks = seq(1, 8))+
  scale_y_continuous(breaks = seq(0, 10, 2))+
  labs(x='\\# species', 
       y='\\# observations')

tikz(file = 'plot_nspeciesDist.tex', standAlone = FALSE, width = figure.width, height = figure.height / 2)
nspeciesDist
dev.off()

# In one subject, only 1 species was found. None of the subjects tested positive on all 8 species.
# In the following table, the odds of finding a one bacteria species when another has been found were investigated.
# In these case, we treat the detection as a binary variable: abuncance bacteria X > 0.0 or not.
bacteria <- armpit %>%
  select(Corynebacterium.1:Corynebacterium.4, Staphylococcus.1:Staphylococcus.4) 

bacteria_table <- matrix(nrow = 2, ncol = 2)
Fisher_exact_p <- matrix(nrow = 8, ncol = 8)
Fisher_exact_or <- matrix(nrow = 8, ncol = 8)

for(i in 1:8){
  for(j in 1:8){
    FT <- fisher.test(x = factor(ifelse(bacteria[i]>0,1,0), levels = c(1,0)), 
                      y =factor(ifelse(bacteria[j]>0,1,0), levels = c(1,0)))
    Fisher_exact_p[i,j] <- FT$p.value
    Fisher_exact_or[i,j] <- FT$estimate
  }
}
  

rownames(Fisher_exact_p) <- colnames(bacteria)
colnames(Fisher_exact_p) <- colnames(bacteria)
rownames(Fisher_exact_or) <- colnames(bacteria)
colnames(Fisher_exact_or) <- colnames(bacteria)

Fisher_exact_or

Fisher_exact_p <- Fisher_exact_p[c(1,2,3,4,6,7,8), c(1,2,3,4,6,7,8)]
Fisher_exact_or <- Fisher_exact_or[c(1,2,3,4,6,7,8), c(1,2,3,4,6,7,8)]

Fisher_exact_p <- Fisher_exact_p %>%
  melt(.) %>% 
  mutate(significant = factor(ifelse(value < 0.05, TRUE, FALSE))) %>%
  rename(p.value = value)

str(Fisher_exact_p)

Fisher_exact_or <- round(Fisher_exact_or,  2) %>% 
  melt(.) %>%
  rename(odds = value) %>%
  mutate(significant = Fisher_exact_p$significant)

# The odds ratios are represented in this table with a color indication of their magnitudes.
# Based on the Fisher Exact result, the significance of this odds ratio can be evaluated.
# The Odds ratios with a p-value < 0.05 are indicated with black text, the others with grey text.
plot_Fisher_exact <- Fisher_exact_or %>%
  ggplot( aes(Var1, Var2)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = odds)) + # background colours are mapped according to the value column
  geom_text(aes(alpha = significant, label = odds)) + # write the values
  scale_alpha_manual(values = c(0.4, 1))+
  scale_fill_gradient2(low = "midnightblue", 
                       mid = "white", 
                       high = "darkred", 
                       midpoint = 3) + # determine the colour
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 14),
        plot.title = element_text(size=16),
        axis.text.y = element_text(size = 14)) + 
  # ggtitle("Fisher-exact: Odds Ratio") +
  theme(legend.position = 'none') + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="Fisher exact\nodds")


tikz(file = 'plot_Fisher_exact.tex', standAlone = FALSE, width = 8, height = 8)
plot_Fisher_exact
dev.off()

#The plot below shows the Kendall correlation coefficients
speciescor_p <- matrix(nrow = 8, ncol = 8)
speciescor_tau <- matrix(nrow = 8, ncol = 8)
for(i in 1:8){
  for(j in 1:8){
    test <- cor.test(bacteria[,i], bacteria[, j], 
                     method = 'kendall', 
                     conf.level = 0.95, 
                     alternative = 'two.sided')
    speciescor_p[i,j] <- test$p.value
    speciescor_tau[i,j] <- test$estimate
  }
}


rownames(speciescor_p) <- colnames(speciescor_p) <- colnames(bacteria)
rownames(speciescor_tau) <- colnames(speciescor_tau) <- colnames(bacteria)

speciescor_p <- speciescor_p %>% 
  melt(.) %>% 
  mutate(significant = factor(ifelse(value < 0.05, TRUE, FALSE))) %>%
  rename(p.value = value)

str(speciescor_p)

speciescor_tau <- round(speciescor_tau,  2) %>% 
  melt(.) %>%
  rename(corr.coeff = value) %>%
  mutate(significant = speciescor_p$significant)


plot_correlation <- speciescor_tau %>%
  ggplot( aes(Var1, Var2)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = corr.coeff)) + # background colours are mapped according to the value column
  geom_text(aes(alpha = significant, label = corr.coeff)) + # write the values
  scale_fill_gradient2(low = "midnightblue", 
                       mid = "white", 
                       high = "darkred", 
                       midpoint = 0) + # determine the colour
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=45, hjust = 1,vjust=1,size = 14),
        plot.title = element_text(size=16),
        axis.text.y = element_text(size = 14)) + 
  ggtitle("Correlation (Kendall)") + 
  theme(legend.position = 'none') + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="correlation coefficients")

tikz(file = 'plot_correlation.tex', standAlone = FALSE, width = 8, height = 8)
  plot_correlation
dev.off()

# For the presentation on December 18th, 2019 we made a poster. 
# The following code exports the png images that were presented in this poster.
ggsave(file="plot_Boxplot_species.png", plot=Boxplot_bacteria, width=10, height=5)
ggsave(file="plot_correlation.png", plot=plot_correlation, width=5.5, height=5.5)
ggsave(file="plot_Fisher_exact.png", plot=plot_Fisher_exact, width=5.5, height=5.5)
ggsave(file="plot_AgeDistribution.png", plot=AgeDist, width=8, height=5)
ggsave(file="plot_Age.png", plot=plot_Age, width=8, height=5)