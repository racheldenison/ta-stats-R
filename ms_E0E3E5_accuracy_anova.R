# ms_E0E3E5_accuracy_anova.R

# clears workspace:  
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
path="/Local/Users/denison/Google Drive/NYU/Projects/Temporal_Attention/Code/Microsaccade/models"
setwd("/Local/Users/denison/Google Drive/NYU/Projects/Temporal_Attention/Code/Microsaccade/data/combined")

#### load data ####
analysisName = "E0E3E5_msBehavAccuracyNormalized_N30"
data <- read.table(paste(analysisName,".txt",sep=""), header = T, sep=" ") 

data$subject <- factor(data$subject)
data$experiment <- factor(data$experiment)
data$target <- factor(data$target)
data$validity <- factor(data$validity)

# visualize data
plot(data$validity, data$acc)
interaction.plot(data$validity, data$target, data$acc)

#change type of each variable
#dataSet=data.frame(mapply(function(dataSet, class) match.fun(paste("as", class, sep="."))(dataSet), dataSet, colClasses("fffd"), SIMPLIFY=FALSE), stringsAsFactors=FALSE)

#set treatment contrasts
#contrasts(data$cue) <- contr.treatment(levels(data$cue), base=which(levels(data$cue) == "t1"))

lmm <- lmer(acc ~ target * validity * experiment + (1|subject), data=data)
Anova(lmm)
save(lmm, file = paste(path, '/', paste(analysisName,'model',sep='_'), ".rda", sep=""))

# by target
lmmResults <- ddply(data, 'target', function(x){
  ## LMM ##
  lmm <- lmer(acc ~ validity * experiment + (1|subject), data=x)
  print(summary(lmm))
  #save model as R files
  save(lmm, file = paste(path, '/', paste(analysisName,x$target[1],'model',sep='_'), ".rda", sep=""))
  #store resulting statistics and p value in a data frame
  
  ###### STOPPED WORKING HERE #######
  data.frame(beta_validity = summary(lmm)$coefficients[2,1], pval_validity = Anova(lmm)[1,3])
})


#save model results as R files
save(lmmResults, file = paste(path, '/', paste(analysisName,'model',sep='_'), ".rda", sep=""))

bootstrapping = TRUE

binStarts <- unique(data$bin_start)
targets <- unique(data$target)

#combine model names
listOfBootstrappedModels <- c()
for (j in 1:length(targets)){
for (i in 1:length(binStarts)){
  listOfBootstrappedModels <- c(listOfBootstrappedModels, paste(analysisName,targets[j],binStarts[i],'model',sep='_'))
}}

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


