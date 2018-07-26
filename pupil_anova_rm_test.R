# anova_test.R

# compare with TAstats.m

# install.packages("ez")

library(ez)

data <- read.csv("/Volumes/purplab/EXPERIMENTS/1_Current Experiments/Rachel/Temporal_Attention/eyeanalysis/xls/anova_test_E0E3_time86.csv")

data$condition <- factor(data$condition)
data$subject <- factor(data$subject)

anova = ezANOVA(
	data=data, 
	dv=pupil, 
	wid=subject, 
	within=.(condition)
)

# example of multi-factor within-subjects ANOVA
# acc_anova = ezANOVA(
#	data=data, 
#	dv=accuracy, 
#	wid=subject, 
#	within=.(validity,early.late,contrast)
#)