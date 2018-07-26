# clears workspace:  
rm(list=ls(all=TRUE)) 

#function to detach all packages (to avoid interferences)
detachAllPackages <- function() { 
  #store basic packages names in a list
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")  
  #make list of all loaded packages
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  #remove basic packages from the list
  package.list <- setdiff(package.list,basic.packages)
  #remove all packages from the list
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
#detach all packages
#detachAllPackages()


###### general settings #######
#save path to main directory
#path="/Users/rachel/Google Drive/NYU/Projects/Temporal_Attention/Code/R/LMM/data32"
path="/Local/Users/denison/Google Drive/NYU/Projects/Temporal_Attention/Code/R/LMM/data32"

# analysis settings
analysisName <- "beta_t1t2vni_cf"
contrastType <- "sum" # "treatment","sum"
comparisonCond <- "valid" # only applies if contrastType is treatment
columnClasses <- "fffffd" #"fffd", "ffffd", "fffffd" make sure ezANOVA and LMM are correct
attentionNames <- c("valid",    "neutral",    "invalid") #c("attended", "neutral", "unattended") #c("valid",    "neutral",    "invalid")

#bootstrap statistical models
bootstrapping <- TRUE


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

#set general theme to black and white
theme_set(theme_bw(10))
#vector with colours
myColours <- c("#7A7A7A","#CC2929", "#61A220", '#0F6499')

###### data collection and transformation #######
print("####################### data collection and transformation #############################")

#change working directory to data directory
setwd(path)
#read in data from txt file
dataSet.wide <- read.table(paste(analysisName,".txt",sep=""), header = T, sep="\t") 

###### split condition variables #######
print("####################### split condition variables #############################")

#add subject name to dataSet, extract from combined variable
dataSet.wide$subNo <- substr(dataSet.wide$subNoxExperiment, start = 1, stop = 2)

#add experiment name to dataSet, extract from combined variable
dataSet.wide$experiment <- substr(dataSet.wide$subNoxExperiment, start = 3, stop = 4)

#delete column
dataSet.wide <- subset(dataSet.wide, select = - subNoxExperiment)

###### reshape dataset #######
#change data from wide to long format
dataSet <- melt(dataSet.wide, measure.vars = attentionNames, value.name = 'pupilMeasure', variable.name = "validity")


###### data preparation #######
print("####################### data preparation #############################")

#change type of each variable
dataSet=data.frame(mapply(function(dataSet, class) match.fun(paste("as", class, sep="."))(dataSet), dataSet, colClasses(columnClasses), SIMPLIFY=FALSE), stringsAsFactors=FALSE)


##### Analysis of latencies #####
print("####################### Analysis of latenciess #############################")

#traditional ANOVA
#ezANOVA(data = dataSet, dv = pupilMeasure, wid = subNo, within = .(validity), between=experiment)
#ezANOVA(data = dataSet, dv = pupilMeasure, wid = subNo, within = .(validity, target), between=experiment)
#ezANOVA(data = dataSet, dv = pupilMeasure, wid = subNo, within = .(validity, target, correctness), between=experiment)

if (contrastType=="sum"){
  # sum contrasts for all
  options(contrasts=c('contr.sum','contr.poly'))
  contrastStr = "sum"
} else if (contrastType=="treatment") {
  #set treatment contrasts
  contrasts(dataSet$validity) <- contr.treatment(levels(dataSet$validity),
                                                 base=which(levels(dataSet$validity) == comparisonCond))
  contrastStr = comparisonCond
}

## LMM ##
#print(summary(pupilMeasure.lmm <- lmer(pupilMeasure ~ validity * experiment + (1|subNo),
#                                                   data=dataSet)))
#print(summary(pupilMeasure.lmm <- lmer(pupilMeasure ~ validity * target * experiment + (1|subNo),
#                                       data=dataSet)))
print(summary(pupilMeasure.lmm <- lmer(pupilMeasure ~ validity * target * correctness * experiment + (1|subNo),
                                       data=dataSet)))

#save model as R files
save(pupilMeasure.lmm, file = paste(path, "/models/", paste(analysisName,contrastStr,sep='_'), ".rda", sep=""))

#do ANOVA equivalent on LMM
Anova(pupilMeasure.lmm)


#combine model names
listOfBootstrappedModels <- c(paste(analysisName,contrastStr,sep='_'))
#listOfBootstrappedModels <- c("latency_t1vni_correct_neutral","latency_t2vni_correct_neutral","beta_t1vni_correct_neutral","beta_t2vni_correct_neutral")
#listOfBootstrappedModels <- c("latency_t1anu_attended","latency_t2anu_attended","beta_t2anu_attended","beta_t2anu_attended")

#only if models should be bootstrapped
if (bootstrapping == TRUE){
  #set number of bootstraps
  n.bootstraps <- 1000
  #set number of bootstraps
  n.bootstrap.reps <- 10
  #loop through the to be bootstrapped models
  for (j in 1:length(listOfBootstrappedModels)){
    #load the respective model
    modelName <- load(paste(path, "/models/", listOfBootstrappedModels[j], ".rda", sep=""))
    #loop through bootstrap circles
    for (i in 1:n.bootstrap.reps){
      #parametric bootstrapping of the the fixed effects
      model.boot <- bootMer(get(modelName), function(x) c(fixef(x)), nsim=n.bootstraps, parallel = "multicore", ncpus = detectCores())
      #read out simulated estimates
      model.boot.estimates <- as.data.frame(model.boot)
      #store file name
      print(filename <- paste(paste(paste(paste(path, "/models/bootstraps/", listOfBootstrappedModels[j], sep=""),"estimates", sep="_"), i*n.bootstraps, sep="_"), ".txt", sep="" ))
      #write estimates as txt file
      write.table(x=model.boot.estimates, file=filename, sep = "\t", row.names = F)
    }
  }

  #change working directory to data directory
  setwd(paste(path, "/models/bootstraps/", sep=""))
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
    print(filename <- paste(paste(paste(path, "/models/", listOfBootstrappedModels[j], sep=""),"BootstrappingResults", sep="_"), ".txt", sep="" ))
    #write estimates as txt file
    write.table(x=bootstrapResults, file=filename, sep = "\t", row.names = F)
  }
}#if (bootstrapping == TRUE){

