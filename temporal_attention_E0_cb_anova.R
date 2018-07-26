# temporal_attention_anova.R

data <- read.csv("TemporalAttention_E0_cb_AllData.csv")

library(ez)

data$subject <- factor(data$subject)
data$contrast <- factor(data$contrast)
data$early.late <- factor(data$early.late)
data$validity <- factor(data$validity)

# visualize data
interaction.plot(data$validity, data$early.late, data$accuracy)

# all data
data0 <- data # store a safe copy of data

acc_anova = ezANOVA(
	data=data, 
	dv=accuracy, 
	wid=subject, 
	within=.(validity,early.late,contrast)
)

# no MA
data_noMA <- subset(data, subject!=7)

# contrasts separately
data_contrast64 <- subset(data, contrast==64)

acc_anova_contrast64 = ezANOVA(
	data=data_contrast64, 
	dv=accuracy, 
	wid=subject, 
	within=.(validity,early.late)
)

data_contrast100 <- subset(data, contrast==100)

acc_anova_contrast100 = ezANOVA(
	data=data_contrast100, 
	dv=accuracy, 
	wid=subject, 
	within=.(validity,early.late)
)

# average across contrasts
data_contrast_ave <- data_contrast64
data_contrast_ave$accuracy <- (data_contrast64$accuracy + data_contrast100$accuracy)/2
data_contrast_ave$RT <- (data_contrast64$RT + data_contrast100$RT)/2

acc_anova_contrast_ave = ezANOVA(
	data=data_contrast_ave, 
	dv=accuracy, 
	wid=subject, 
	within=.(validity,early.late)
)

acc_anova_T1_ave = ezANOVA(
	data=subset(data_contrast_ave, early.late=="early"),
	dv=accuracy, 
	wid=subject, 
	within=.(validity)
)

acc_anova_T2_ave = ezANOVA(
	data=subset(data_contrast_ave, early.late=="late"),
	dv=accuracy, 
	wid=subject, 
	within=.(validity)
)

# early and late separately
data_T1 <- subset(data, early.late=="early")

acc_anova_T1 = ezANOVA(
	data=data_T1, 
	dv=accuracy, 
	wid=subject, 
	within=.(validity,contrast)
)

data_T2 <- subset(data, early.late=="late")

acc_anova_T2 = ezANOVA(
	data=data_T2, 
	dv=accuracy, 
	wid=subject, 
	within=.(validity,contrast)
)

# average across contrasts and targets
data_contrast_ave_T1 <- subset(data_contrast_ave, early.late=="early")
data_contrast_ave_T2 <- subset(data_contrast_ave, early.late=="late")
data_contrast_target_ave <- data_contrast_ave_T1
data_contrast_target_ave$accuracy <- (data_contrast_ave_T1$accuracy + data_contrast_ave_T2$accuracy)/2
data_contrast_target_ave$RT <- (data_contrast_ave_T1$RT + data_contrast_ave_T2$RT)/2

# t-tests
pairwise.t.test(data_contrast_target_ave$accuracy, data_contrast_target_ave$validity, p.adj="holm", paired=TRUE)

pairwise.t.test(data_contrast_ave_T1$accuracy, data_contrast_ave_T1$validity, p.adj="holm", paired=TRUE)

pairwise.t.test(data_contrast_ave_T2$accuracy, data_contrast_ave_T2$validity, p.adj="holm", paired=TRUE)

### RT
rt_anova = ezANOVA(
	data=data, 
	dv=RT, 
	wid=subject, 
	within=.(validity,early.late,contrast)
)

# t-tests
pairwise.t.test(data_contrast_target_ave$RT, data_contrast_target_ave$validity, p.adj="holm", paired=TRUE)

pairwise.t.test(data_contrast_ave_T1$RT, data_contrast_ave_T1$validity, p.adj="holm", paired=TRUE)

pairwise.t.test(data_contrast_ave_T2$RT, data_contrast_ave_T2$validity, p.adj="holm", paired=TRUE)
