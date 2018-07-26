# ms_E0E3E5_rebound_anova.R

rm(list=ls(all=TRUE)) 

#### load libraries ####
#needed for ezANOVA()
library(ez)
#load R.utils
library(R.utils)
#need for ddply
library(plyr)
#needs this library
library(gdata)
#needs this library
library(ggplot2)
#needed for recode
library(car)
#needed for melt
library(reshape2)
#needed for unit
library(gridExtra)
#needed for gtable
library(gtable)
#needed for glmer
library(lme4)
#needed for string comparison
library(stringr)
#needed for testInteractions
library(phia)
#needed for detectCores()
library(parallel)
#load tidyr
library(tidyr)

#### set directory ####
home="/Local/Users/denison/Google Drive/NYU/Projects"
#home="~/Documents/NYU/Projects"
path=paste(home, "Temporal_Attention/Code/Microsaccade/models", sep="/")
setwd(paste(home, "Temporal_Attention/Code/Microsaccade/data/combined", sep="/"))
#setwd("~/Documents/NYU/Projects/Temporal_Attention/XLS")
#setwd("/Local/Users/denison/Google Drive/NYU/Projects/Temporal_Attention/XLS")

#### load data ####
#analysisName = "rebound"
#data_wide <- read.csv("MS_E0E3E5_Rebound_n30.csv")
analysisName = "E0E3E5_msRebound2Latency_N30_dataZ" #E5_msNearTarget_accNT_-100to0_N9"
data_wide <- read.table(paste(analysisName,".txt",sep=""), header = T, sep=" ")

#change data from wide to long format
data <- melt(data_wide, measure.vars = c("neutral", "t1", "t2"), value.name = 'latency', variable.name = "cue")
#data <- melt(data_wide, measure.vars = c("neutral", "t1", "t2","t3"), value.name = 'latency', variable.name = "cue") # 3 targets

data$subject <- factor(data$subject)
data$experiment <- factor(data$experiment)
data$cue <- factor(data$cue)

### VBehav only ###
data$target <- factor(data$target)
data$accuracy <- factor(data$accuracy)
###################
### msNearTarget ###
data$target <- factor(data$target)
data$ms_in_window <- factor(data$ms_in_window)
###################

# visualize data
plot(data$cue, data$latency)

#change type of each variable
#dataSet=data.frame(mapply(function(dataSet, class) match.fun(paste("as", class, sep="."))(dataSet), dataSet, colClasses("fffd"), SIMPLIFY=FALSE), stringsAsFactors=FALSE)

#set treatment contrasts
contrasts(data$cue) <- contr.treatment(levels(data$cue), base=which(levels(data$cue) == "t2")) # t1

## LMM ##
print(summary(lmm <- lmer(latency ~ cue * experiment + (1|subject), data=data)))
#print(summary(lmm <- lmer(latency ~ cue + (1|subject), data=data))) # single experiment
print(summary(lmm <- lmer(latency ~ cue * experiment * target * accuracy + (1|subject), data=data))) # vs. behav

x = subset(data, target=="t3")
print(summary(lmm <- lmer(acc ~ ms_in_window * experiment + (1|subject), data=x))) # ms near target - accNT
print(summary(lmm <- lmer(ms_prop ~ cue * experiment + (1|subject), data=x))) # ms near target - msNT

#save model as R files
save(lmm, file = paste(path, '/', paste(analysisName,'baset2_model',sep='_'), ".rda", sep=""))

#do ANOVA equivalent on LMM
Anova(lmm)

bootstrapping = TRUE

#combine model names
listOfBootstrappedModels <- c(paste(analysisName,'baset2_model',sep='_'))

#only if models should be bootstrapped
if (bootstrapping == TRUE){
  #set number of bootstraps
  n.bootstraps <- 1000
  #set number of bootstraps
  n.bootstrap.reps <- 10
  #loop through the to be bootstrapped models
  for (j in 1:length(listOfBootstrappedModels)){
    #load the respective model
    modelName <- load(paste(path, "/", listOfBootstrappedModels[j], ".rda", sep=""))
    #loop through bootstrap circles
    for (i in 1:n.bootstrap.reps){
      #parametric bootstrapping of the the fixed effects
      model.boot <- bootMer(get(modelName), function(x) c(fixef(x)), nsim=n.bootstraps, parallel = "multicore", ncpus = detectCores())
      #read out simulated estimates
      model.boot.estimates <- as.data.frame(model.boot)
      #store file name
      print(filename <- paste(paste(paste(paste(path, "/bootstraps/", listOfBootstrappedModels[j], sep=""),"estimates", sep="_"), i*n.bootstraps, sep="_"), ".txt", sep="" ))
      #write estimates as txt file
      write.table(x=model.boot.estimates, file=filename, sep = "\t", row.names = F)
    }
  }
  
  #change working directory to data directory
  setwd(paste(path, "/bootstraps/", sep=""))
  #loop through the to be bootstrapped models
  for (j in 1:length(listOfBootstrappedModels)){
    
    #list all txt files in data directory
    listOfFiles.boot = list.files(pattern= paste(listOfBootstrappedModels[j], "_estimates",sep=""))
    #read every .txt file in and append it to data frame
    bootstraps  <-  do.call(rbind, lapply(listOfFiles.boot, function(fname) {
      #print to console the name of the file
      print(fname)
      #read in the file
      dum  <-  read.table(fname,header = T)
      #return data.frame to calling function
      return(dum)
    }))
    #reshape bootstraps to long format
    bootstraps.long <- gather(bootstraps, fixedEffect, estimate)
    
    #calculate confidence intervals from quantiles
    bootstrapResults <- ddply(bootstraps.long, "fixedEffect", function(x) {
      #read out quantiles for CI
      tmp <- quantile(x$estimate, probs=c(0.05,0.95), na.rm=T)
      #proportion of sims in which parameter estimate was larger than zero
      prop <- sum(x$estimate>0, na.rm=T)/sum(!is.na(x$estimate))
      # compare simulated proportion against chance proportion (from pvals.fnc)
      pValue <- 2 * pmax(0.5/sum(!is.na(x$estimate)), pmin(prop, 1 - prop))
      #combine everything into one data frame
      data.frame(meanEstimate=round(mean(x$estimate, na.rm=T),2), CI.low=round(tmp[1],2), CI.high=round(tmp[2],2), prop=round(prop,3), pValue=round(pValue,3))
    })
    
    #store file name
    print(filename <- paste(paste(paste(path, "/", listOfBootstrappedModels[j], sep=""),"BootstrappingResults", sep="_"), ".txt", sep="" ))
    #write estimates as txt file
    write.table(x=bootstrapResults, file=filename, sep = "\t", row.names = F)
  }
}#if (bootstrapping == TRUE){


