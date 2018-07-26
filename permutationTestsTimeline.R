# clears workspace:  
rm(list=ls(all=TRUE)) 

#### cluster permutation tests #####
print('#### cluster permutation tests #####')

#load R.utils
library(R.utils)
#needed for recode
library(car)
#load plyr for ddply()
library(plyr)
#load reshape2 for melt()
library(reshape2)
#needed for unit
library(grid)
#needed for lmer
library(lme4)
#for parallel use of cores
library(doMC)
#register all four cores
registerDoMC(3)
#parallel for loop
library(foreach)
# graphics
library(ggplot2)

#home="~/Documents/NYU/Projects"
home="/Local/Users/denison/Google Drive/NYU/Projects"
path=paste(home,"Temporal_Attention/Code/Microsaccade/data/combined/",sep="/")

# analysis settings
analysisName <- "E0E3E5_msrate_-500to0ms_N30"
attentionNames <- c("neutral","t1","t2")
contrastType <- "treatment" # "treatment","sum"
comparisonCond <- "t1" # only applies if contrastType is treatment
columnClasses <- "ffdfd" # make sure LMM is correct
measureName <- "msrate"
betaName <- "beta_cueneutral"
nSubjects <- 30

#savefilename = paste('TimelinesPermutationTests', measureName, sep='_')
savefilename = paste('TimelinesPermutationTests', analysisName, sep='_')


#### get data ####
#change working directory to data directory
setwd(path)
#read in data from txt file
dataSet.wide <- read.table(paste(analysisName,".txt",sep=""), header = T, sep=" ") 

#change data from wide to long format
dataSet <- melt(dataSet.wide, measure.vars = attentionNames, value.name = measureName, variable.name = "cue")

#change type of each variable
dataSet=data.frame(mapply(function(dataSet, class) match.fun(paste("as", class, sep="."))(dataSet), dataSet, colClasses(columnClasses), SIMPLIFY=FALSE), stringsAsFactors=FALSE)

if (contrastType=="sum"){
  # sum contrasts for all
  options(contrasts=c('contr.sum','contr.poly'))
  contrastStr = "sum"
} else if (contrastType=="treatment") {
  #set treatment contrasts
  contrasts(dataSet$cue) <- contr.treatment(levels(dataSet$cue),
                                            base=which(levels(dataSet$cue) == comparisonCond))
  contrastStr = comparisonCond
}


#function to compare conditions using paired t-test for each timepoint
#as this is a within design, the data are not permuted freely, 
#but rather for each subject the two conditions are randomly assigned
#to do so, the data.frame is sorted by condition and then by subject
#then datapoints are assigned to the two conditions based on variable subSet, that is supposed to be a vector of
#zeros and ones exactly half as long as the data frame
timepointwisePermutationTests <- function(dataSet, cond){
  #divide input data by timepoint
  ddply(dataSet, 'time', function(x){
    # print(x$time[1])
    # order x experiemnt, subject, then cue
    x <- x[with(x, order(experiment, subject, cue)), ]
    # attach condition column - might be permuted
    x$cond <- factor(cond)
    contrasts(x$cond) <- contr.treatment(levels(x$cond),
                                              base=which(levels(x$cond) == "2"))
    # print(x)
    # print("-----")
    ## LMM ##
    lmm <- lmer(msrate ~ cond * experiment + (1|subject), data=x)
    #store resulting statistic and p value in a data frame
    coeff <- summary(lmm)$coefficients
    data.frame(Int = coeff[1,1], beta_cueneutral = coeff[2,1], beta_cuet2 = coeff[3,1], beta_e3 = coeff[4,1], beta_e5 = coeff[5,1],
               beta_cueneutral_e3 = coeff[6,1], beta_cuet2_e3 = coeff[7,1], beta_cueneutral_e5 = coeff[8,1],
               beta_cuet2_e5 = coeff[9,1], pval_cue = Anova(lmm)[1,3])
    #data.frame(t.Value = ttestResults$statistic, p.Value = ttestResults$p.value)
  }, .parallel = F)
}

# nPermutations <- 20 #1000
if (!file.exists(paste(path, savefilename, '.txt', sep=''))){
  #start by running the t-tests on the original data
  dataStat <- cbind(timepointwisePermutationTests(dataSet,
                  foreach(i=1:nSubjects, .combine=rbind) %do% {as.matrix(1:3)}), permutationNo = 0)
  
  #now run on permuted data
  for (j in 160:1000){
    print(j)
    cond <- foreach(i=1:nSubjects,.combine=rbind) %do% {as.matrix(sample(1:3, 3))}
    write.table(cond, paste(path, 'cond.txt', sep=''))
    newStat <- cbind(timepointwisePermutationTests(dataSet, cond), permutationNo = j)
    dataStat <- rbind(dataStat, newStat)
    write.table(dataStat, paste(path, savefilename, '.txt', sep=''))
    if (j %% 10==0){
      print("sleep")
      Sys.sleep(10)
    }
  }

  # dataStat <- rbind(dataStat, foreach(i=1:nPermutations,.combine=rbind) %dopar% {cbind(timepointwisePermutationTests(dataSet,
  #            foreach(i=1:nSubjects,.combine=rbind) %dopar% {as.matrix(sample(1:3, 3))}), permutationNo = i)})
  # #save test results of permuted data
  # write.table(dataStat, paste(path, savefilename, '.txt', sep=''))
}

#read in the permutations
dataStat <- read.table(paste(path, savefilename, '.txt', sep=''))
nPermutations <- length(unique(dataStat$permutationNo))-1

#check for significance and identify clusters
#find 5% percentile for each timepoint (but exclude original data)
tmp <- ddply(dataStat, 'time', function(x){
  print(x$time[1])
  x$significant.permtested.neg <- as.numeric(x[,betaName] < quantile(x[x$permutationNo != 0,betaName], 0.025)) 
  x$significant.permtested.pos <- as.numeric(x[,betaName] > quantile(x[x$permutationNo != 0,betaName], 0.975))
  #x$significant.permtested.neg <- as.numeric(x$beta_cueneutral < quantile(x$beta_cueneutral[x$permutationNo != 0], 0.025)) 
  #x$significant.permtested.pos <- as.numeric(x$beta_cueneutral > quantile(x$beta_cueneutral[x$permutationNo != 0], 0.975))
  x#return x
})
#merge x with dataStat
relevantCols <- c('time', 'permutationNo', betaName)
dataStatBeta <- merge(dataStat[,relevantCols], tmp[,c(relevantCols, "significant.permtested.pos", "significant.permtested.neg")], by = relevantCols)
#delete tmp
rm(tmp)

# check for significance for each beta value
# find 5% percentile for each timepoint (but exclude original data)
# betaNames = c("beta_cueneutral", "beta_cuet2", "beta_e3", "beta_e5", "beta_cueneutral_e3",
#                "beta_cuet2_e3", "beta_cueneutral_e5", "beta_cuet2_e5")
# for (betaName in betaNames){
#   print(betaName)
# dataStat <- ddply(dataStat, 'time', function(x){
#   x[,paste(betaName, ".significant.permtested.neg", sep="")] <- as.numeric(x[,betaName] < quantile(x[x$permutationNo != 0,betaName], 0.025))
#   x[,paste(betaName, ".significant.permtested.pos", sep="")] <- as.numeric(x[,betaName] > quantile(x[x$permutationNo != 0,betaName], 0.975))
#   #x$significant.permtested.neg <- as.numeric(x$beta_cueneutral < quantile(x$beta_cueneutral[x$permutationNo != 0], 0.025))
#   #x$significant.permtested.pos <- as.numeric(x$beta_cueneutral > quantile(x$beta_cueneutral[x$permutationNo != 0], 0.975))
#   x#return x
# })}

#determine clusters for each permutation seperately
clusterStat <- ddply(dataStatBeta, 'permutationNo', function(x){
  #sort by time
  x <- x[order(x$time),]
  #label significant timepoints with the consecutive cluster number
  x$clusterNo.pos <- cumsum(diff(c(0, x$significant.permtested.pos)) == 1) * x$significant.permtested.pos
  #label significant timepoints with the consecutive cluster number
  x$clusterNo.neg <- cumsum(diff(c(0, x$significant.permtested.neg)) == 1) * x$significant.permtested.neg
  #merge positive and negative clusters
  x$clusterNo <- 1000 * x$clusterNo.pos + x$clusterNo.neg
  #split timeline into clusters
  ddply(subset(x, clusterNo > 0), 'clusterNo', function(z){
    #export sum of t values, start point, end point
    z$clusterStat.Summed <- sum(z[,betaName])
    z$timepoint.start <- min(z$time)
    z$timepoint.end <- max(z$time)
    z
  })
}, .parallel = F)

#select cluster with highest sum of t-values
clusterStatSum <- ddply(clusterStat, 'permutationNo', function(x){
  #choose only the cluster with maximal summed t value
  x[which.max(abs(x$clusterStat.Summed)),]
}, .parallel = F)

#expand grid as combinations without any significant comparison are missing
clusterStatSum <- merge(clusterStatSum, data.frame(permutationNo = unique(dataStat[,'permutationNo'])), 
              by='permutationNo', all = T)
#replace NA with 0
clusterStatSum$clusterStat.Summed[is.na(clusterStatSum$clusterStat.Summed)] <- 0

#select cluster with highest sum of t-values
#compare sum of t values from original set with permuted sums
print(data.frame(p.Value.ClusterSize = mean(clusterStatSum$clusterStat.Summed[clusterStatSum$permutationNo != 0] > clusterStatSum$clusterStat.Summed[clusterStatSum$permutationNo == 0]),
           timepoint.start = clusterStatSum$timepoint.start[clusterStatSum$permutationNo == 0], timepoint.end = clusterStatSum$timepoint.end[clusterStatSum$permutationNo == 0]))

# find 95% interval for max cluster sum
threshNeg <- quantile(clusterStatSum$clusterStat.Summed[1:nPermutations+1], 0.025) # should this exclude permutation==0?
threshPos <- quantile(clusterStatSum$clusterStat.Summed[1:nPermutations+1], 0.975)

# find clusters that exceed the threshold values
x = subset(clusterStat, permutationNo==0)
x = unique(x[,(ncol(x)-2):ncol(x)])
sigClusters.pos <- x[x$clusterStat.Summed > threshPos,]
sigClusters.neg <- x[x$clusterStat.Summed < threshNeg,]

sigClusters <- merge(sigClusters.neg, sigClusters.pos, all=TRUE, sort=FALSE)

# calculate p-values
for (i in 1:nrow(sigClusters)){
  sigClusters$pval[i] = mean(clusterStatSum$clusterStat.Summed[clusterStatSum$permutationNo != 0] > sigClusters$clusterStat.Summed[i])
}

print(sigClusters)

write.table(sigClusters, paste(path, savefilename, '_', betaName, '.txt', sep=''))

# plot
# gp <- ggplot(dataStat, aes(x=time, y=beta_cuet2, colour = permutationNo))
# gp <- gp + geom_point()
# gp <- gp + geom_hline(aes(yintercept = 0))
# gp

# save perm0 (real data)
#dataStat0 = subset(dataStat, permutationNo==0)
#dataStat0 <- dataStat0[order(dataStat0$time),]
#write.table(dataStat0, paste(path, savefilename, '_perm0.txt', sep=''))

# check number of elements
# a = c()
# for (i in 1:1000){a[i] = nrow(subset(dataStat, permutationNo==i))}
