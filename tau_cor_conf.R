#################################################################################################
# Project Principles of Statistical data analysis                                               #
# 2019 - 2020                                                                                   #
# Group 9 : Jan Alexander, Paul Morbee, Joren Verbeke, Steven Wallaert                          #
#                                                                                               #
# Function tau_cor_conf                                                                         #
# Version: 1.0      Date: 16/12/2019                                                            #
#                                                                                               #
# Functionality:    Calculates Kendall's tau statistic, p value and confidence interval         #
#                   Based on the R corr function                                                #
#                   P value based on a z-test described in literature (sample size at least 10) #
#                   Confidence interval calculated with bootstrapping                           #
#                                                                                               #
# Input named vars: df:data frame containing the columns to be associated                       #
#                   C1: Column name or number of first variable                                 #
#                   C2: Column name or number of second variable                                #
#                   repl: number of simulations (default=2000)                                  #
#                   conf: desired confidence level (default=0.95)                               #
#                                                                                               #
# Output: a list with following components:                                                     #
#                   $ first_column : chr    First column of the correlaion                      #
#                   $ second_column: chr    Second column of the correlaion                     #
#                   $ tau.value    : num    Calculated tau value                                #  
#                   $ tau.pvalue   : num    pvalue for H0 tau=0 versus Ha tau<>0                #
#                   $ conf_level   : num    Chosen confidence level (default=0.95)              #
#                   $ conf_lower   : num    Confidence Interval (lower)                         #
#                   $ conf_upper   : num    Confidence Interval (upper)                         #
#                                                                                               #
# Errors:           Erroneous data input                                                        #
# Warnings:         Whenever the CI cannot be calculated (equal colums)                         #
#################################################################################################


library(boot)
tau_cor_conf <- function(df,C1,C2, repl=2000, conf=0.95){

# Control of the df parameter
  if(class(df)!="data.frame") stop("df parameter is not a data frame")

# Control of the C1 parameter
  if(length(C1)!=1) stop("One and only one column allowed in C1")

# C1 can be either the number of a column or its name  
# C1 must be a valid name or column number  
  if(is.numeric(C1)) {            
    if(C1<1 | C1>ncol(df)) {
      stop("C1 parameter is not the number of a column")
    } else {
    num_col1 <- C1
    } 
  } else {
    if(is.character(C1)) {
     if(sum(names(df)==C1)!=1) {
       stop("C1 is not a valid column name") 
     } else {
       num_col1 <- which(names(df)==C1)
     }
    } else {
      stop("C1 parameter is neither a column number nor a column name")
    }
  }

# Control of the C2 parameter  

  if(length(C2)!=1) stop("One and only one column allowed in C2")

# C2 can be either the name of a column or its name  
# C2 must be a valid name or column number  
  
  if(is.numeric(C2)) {
    if(C2<1 | C2>ncol(df)) {
      stop("C2 parameter is not the number of a column")
    } else {
      num_col2 <- C2
    } 
  } else {
    if(is.character(C2)) {
      if(sum(names(df)==C2)!=1) {
        stop("C2 is not a valid column name") 
      } else {
        num_col2 <- which(names(df)==C2)
      }
    } else {
      stop("C2 parameter is neither a column number nor a column name")
    }
  }

# Control the repl variable  
  if(length(repl)!=1) stop("One and only one value in repl")  
  if(!is.numeric(repl)) stop("repl parameter is not numeric")
  repl=abs(as.integer(repl))
  if(repl==0) repl <- 2000

# Control the conf variable  
  if(length(conf)!=1) stop("One and only one value in conf")  
  if(!is.numeric(conf)) stop("conf parameter is not numeric")
  if(conf<0 | conf>1) stop("conf parameter must be between 0 and 1")

# C1 and C2 must be ordered data so no character vector nor a not ordered factor  
  if(is.character(df[,num_col1])) stop("C1 is character vector, must be an ordered factor")
  if(is.factor(df[,num_col1]) & !is.ordered(df[,num_col1])) stop("Factor in C1 must be ordered")
  
  if(is.character(df[,num_col2])) stop("C2 is character vector, must be an ordered factor")
  if(is.factor(df[,num_col2]) & !is.ordered(df[,num_col2])) stop("Factor in C2 must be ordered")
  
#first of all: bootstrapping doesn't handle equal vectors so we have to treat this case separately
  if(identical(df[,num_col1],df[,num_col2])) {
# In case of equal vectors we already fill in the list object (for instance tau=1)
    warning("Identical columns in C1 and C2, CI can not be calculated with bootstrapping")
    tau.value <- 1
    conf_level <- conf
    conf_lower <- NA
    conf_upper <- NA
  } else {
#We construct the bootstrapping
#We repeat a call to the bootTau function and simulate different corr calculations
#with samples based on the data frame. Boot creates an object (boot_kendall)   
    boot_kendall <- boot(df, bootTau, repl, C1=num_col1, C2=num_col2)

#based on the boot_kendall object we create the boot.ci object   
    obj_boot.ci  <- boot.ci(boot_kendall, type="norm", conf=conf)

#we extract the different components of the output list  
    tau.value <- obj_boot.ci$t0
    conf_level <- obj_boot.ci$normal[1]
    conf_lower <- obj_boot.ci$normal[2]
    conf_upper <- obj_boot.ci$normal[3]
  }
  first_column <- names(df)[num_col1]
  second_column <- names(df)[num_col2]
  n <- nrow(df)
  if (n < 10) {
    tau.pvalue <- NA      #with sample size smaller then 10 we cannot approximate p value with z-test 
  } else {
    z_tau <- abs(tau.value/sqrt(2*(2*n+5))*3*sqrt(n*(n-1)))  #applying the z-test
    tau.pvalue <- 2*(1-pnorm(z_tau))
  }

#Creation of the ouput list  
  list(first_column=first_column, second_column=second_column,
       tau.value=tau.value, tau.pvalue=tau.pvalue,
       conf_level=conf_level,conf_lower=conf_lower,conf_upper=conf_upper)
}

#Base fonction for the Kendall bootstrapping
bootTau<-function(df,i, C1, C2) {
  if(class(df)!="data.frame") stop("df parameter is not a data frame")
  if(!is.numeric(C1)) stop("C1 parameter of function bootTau is not numeric")
  if(!is.numeric(C2)) stop("C2 parameter of function bootTau is not numeric")
  if(C1<1 | C1>ncol(df)) stop("C1 parameter of function bootTau is not the number of a column")
  if(C2<1 | C2>ncol(df)) stop("C2 parameter of function bootTau is not the number of a column")
  cor(as.numeric(df[i,C1]),as.numeric(df[i,C2]), use="complete.obs", method="kendall")
}
