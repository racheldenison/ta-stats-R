# temporal_attention_E3_adjust_anova.R

setwd("/Users/rachel/Documents/NYU/Projects/Temporal_Attention/XLS")
#data <- read.csv("TemporalAttention_E3_adjust_AllData_mixtureNoBiasMAP.csv")
data <- read.csv("TemporalAttention_E3_adjust_AllData_mixtureNoBiasMAP_session1.csv")

library(ez)

data$subject <- factor(data$subject)
data$T1T2 <- factor(data$T1T2)
data$validity <- factor(data$validity)

# visualize data
interaction.plot(data$validity, data$T1T2, data$absMu)
interaction.plot(data$validity, data$T1T2, data$mu)
interaction.plot(data$validity, data$T1T2, data$g)
interaction.plot(data$validity, data$T1T2, data$sd)
interaction.plot(data$validity, data$T1T2, data$B)

# all data
data0 <- data # store a safe copy of data

absmu_anova = ezANOVA(
	data=data, 
	dv=absMu, 
	wid=subject, 
	within=.(validity,T1T2)
)

mu_anova = ezANOVA(
	data=data, 
	dv=mu, 
	wid=subject, 
	within=.(validity,T1T2)
)

g_anova = ezANOVA(
	data=data, 
	dv=g, 
	wid=subject, 
	within=.(validity,T1T2)
)

sd_anova = ezANOVA(
	data=data, 
	dv=sd, 
	wid=subject, 
	within=.(validity,T1T2)
)

B_anova = ezANOVA(
	data=data, 
	dv=B, 
	wid=subject, 
	within=.(validity,T1T2)
)





# T1 and T2 separately
#data_T1 <- subset(data, T1T2=="T1")
data_T1 <- subset(data, T1T2==1)

absmu_anova_T1 = ezANOVA(
	data=data_T1, 
	dv=absMu, 
	wid=subject, 
	within=.(validity)
)

mu_anova_T1 = ezANOVA(
	data=data_T1, 
	dv=mu, 
	wid=subject, 
	within=.(validity)
)

g_anova_T1 = ezANOVA(
	data=data_T1, 
	dv=g, 
	wid=subject, 
	within=.(validity)
)

sd_anova_T1 = ezANOVA(
	data=data_T1, 
	dv=sd, 
	wid=subject, 
	within=.(validity)
)

B_anova_T1 = ezANOVA(
	data=data_T1, 
	dv=B, 
	wid=subject, 
	within=.(validity)
)



#data_T2 <- subset(data, T1T2=="T2")
data_T2 <- subset(data, T1T2==2)

absmu_anova_T2 = ezANOVA(
	data=data_T2, 
	dv=absMu, 
	wid=subject, 
	within=.(validity)
)

mu_anova_T2 = ezANOVA(
	data=data_T2, 
	dv=mu, 
	wid=subject, 
	within=.(validity)
)

g_anova_T2 = ezANOVA(
	data=data_T2, 
	dv=g, 
	wid=subject, 
	within=.(validity)
)

sd_anova_T2 = ezANOVA(
	data=data_T2, 
	dv=sd, 
	wid=subject, 
	within=.(validity)
)

B_anova_T2 = ezANOVA(
	data=data_T2, 
	dv=B, 
	wid=subject, 
	within=.(validity)
)


# Let's try a randomization test
scramble <- function(x, k=3L) {
    x.s <- seq_along(x)
    y.s <- sample(x.s)
    x[unlist(split(x.s[y.s], (y.s-1) %/% k), use.names = FALSE)]
}

# set the original anova F and the dv below to the variable of interest (g, sd, etc)
# set the original anova F to the entire data set or subset of interest
# set datar to be the entire data set or subset of interest
# set locality of scramble
# set within factors
a0F = B_anova_T2$ANOVA$F # original (observed) F
nF = length(a0F)
N <- 1000
F <- matrix(nrow=N, ncol=nF)
datar <- data_T2 # data # data_T1 # data_T2
x <- c(1:nrow(datar))

for (i in 1:N) {
	# shuffle condition labels within subject
	label_order <- scramble(x,3) # 6 for whole data set, 3 for T1/T2 subset
	datar[,2] <- datar[label_order,2]
	datar[,3] <- datar[label_order,3]
	
	#datars <- datar[with(datar, order(subject, T1T2, validity)),]

	# do anova on shuffled data
	aa <- ezANOVA(
		data=datar, 
		dv=B, 
		wid=subject, 
		within=.(validity) # (validity,T1T2) # (validity)
	)
	F[i,] <- aa$ANOVA$F
}

for (i in 1:nF) {
	print(sum(F[,i]>a0F[i])/N)
}

# visualize F distribution generated by randomization
hist(F[,1],20)
plot(density(F[,1]))
