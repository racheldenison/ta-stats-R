# temporal_attention_anova.R

setwd('~/Documents/NYU/Projects/Temporal_Attention/XLS/')

data <- read.csv("TemporalAttention_E5_T3_AllData.csv")

library(ez)
library(plyr)

data$subject <- factor(data$subject)
data$target <- factor(data$target)
data$validity <- factor(data$validity)

# visualize data
interaction.plot(data$validity, data$target, data$accuracy)

# all data
data0 <- data # store a safe copy of data

acc_anova = ezANOVA(
	data=data, 
	dv=accuracy, 
	wid=subject, 
	within=.(validity,target)
)

data_T1=subset(data, target=="T1")
data_T2=subset(data, target=="T2")
data_T3=subset(data, target=="T3")

data_target_ave <- ddply(data,c("subject","validity"),summarize,accuracy=mean(accuracy),RT=mean(RT))

# t-tests
pairwise.t.test(data_target_ave$accuracy, data_target_ave$validity, p.adj="holm", paired=TRUE)

pairwise.t.test(data_T1$accuracy, data_T1$validity, p.adj="holm", paired=TRUE)
pairwise.t.test(data_T2$accuracy, data_T2$validity, p.adj="holm", paired=TRUE)
pairwise.t.test(data_T3$accuracy, data_T3$validity, p.adj="holm", paired=TRUE)

### RT
rt_anova = ezANOVA(
	data=data, 
	dv=RT, 
	wid=subject, 
	within=.(validity,target)
)

# t-tests
pairwise.t.test(data_target_ave$RT, data_target_ave$validity, p.adj="holm", paired=TRUE)

pairwise.t.test(data_T1$RT, data_T1$validity, p.adj="holm", paired=TRUE)
pairwise.t.test(data_T2$RT, data_T2$validity, p.adj="holm", paired=TRUE)
pairwise.t.test(data_T3$RT, data_T3$validity, p.adj="holm", paired=TRUE)
