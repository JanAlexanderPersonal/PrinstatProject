#################################################################################################
# Project Principles of Statistical data analysis                                               #
# 2019 - 2020                                                                                   #
# Group 9 : Jan Alexander, Paul Morbee, Joren Verbeke, Steven Wallaert                          #
#                                                                                               #
# Testscript test_tau_cor_conf                                                                  #
#                                                                                               #
# Test Object:                                                                                  #
# Function tau_cor_conf                                                                         #
# Version: 1.0      Date: 16/12/2019                                                            #
#                                                                                               #
# Functionality:    This testscript builds a named matrix with all cross correlations in        #
#                   the data frame armpit based on tau_cor_conf                                 #
#                   Tested: Numeric, Factor, Ordered, Character input                           #
#                   Confidence levels and p.value tested on individual cases                    #
#                   Tests integration with source of tau_cor_conf function                      #
#                                                                                               #
# Pre condition:    The armpit data frame should be loaded                                      #
#################################################################################################

source("tau_cor_conf.R")
df<-armpit
df$Gender<-factor(df$Gender,ordered = TRUE)
corr.matrix <- matrix(data=NA, ncol = length(df), nrow = length(df))
row.names(corr.matrix) <- names(df)
colnames(corr.matrix) <- names(df)
for(i in 1:ncol(df)) {
  NA_ind <- FALSE
  if ((is.factor(df[,i]) & !is.ordered(df[,i])) | is.character(df[,i])) NA_ind <- TRUE
  for (j in 1:i) {
    if ((is.factor(df[,j]) & !is.ordered(df[,j]))  | is.character(df[,j]) | NA_ind) {
      corr.matrix[i,j] <- NA
    } else {
      corr.matrix[i,j] <- tau_cor_conf(df,i,j)$tau.value
    }
    corr.matrix[j,i] <- corr.matrix[i,j]
  }  
}

