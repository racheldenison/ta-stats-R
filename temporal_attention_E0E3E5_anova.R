# temporal_attention_anova.R

setwd("~/Documents/NYU/Projects/Temporal_Attention/XLS")
#data <- read.csv("TemporalAttention_E0E3E5_zscored.csv")
data <- read.csv("TemporalAttention_E0E5.csv")

library(ez)

data$subject <- factor(data$subject)
data$exp <- factor(data$exp)
data$target <- factor(data$target)
data$validity <- factor(data$validity)

# visualize data
interaction.plot(data$validity, data$target, data$acc)

# anova, z-score
z_anova = ezANOVA(
	data=data, 
	dv=zscore, 
	wid=subject, 
	within=.(validity,target),
	between=.(exp)
)

# anova, acc
acc_anova = ezANOVA(
	data=data, 
	dv=acc, 
	wid=subject, 
	within=.(validity,target),
	between=.(exp)
)

# without E3

# targets separately
data_T1 = subset(data, target==1)
data_T2 = subset(data, target==2)

# t-tests, z-score
pairwise.t.test(data$zscore, data$validity, p.adj="none", paired=TRUE)

pairwise.t.test(data_T1$zscore, data$validity, p.adj="none", paired=TRUE)
pairwise.t.test(data_T2$zscore, data$validity, p.adj="none", paired=TRUE)

# t-tests, acc
pairwise.t.test(data$acc, data$validity, p.adj="none", paired=TRUE)

pairwise.t.test(data_T1$acc, data$validity, p.adj="none", paired=TRUE)
pairwise.t.test(data_T2$acc, data$validity, p.adj="none", paired=TRUE)



