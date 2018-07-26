# temporal_attention_E3_adjust_nonparametric_anova.R

setwd("/Users/rachel/Documents/NYU/Projects/Temporal_Attention/XLS")

# http://depts.washington.edu/aimgroup/proj/art/
library(ARTool)

## using the ARTool library -- gives us partial eta squared
data <- read.csv("TemporalAttention_E5_T3_AllData_rt.csv")

data$subject <- factor(data$subject)
data$target <- factor(data$target)
data$validity <- factor(data$validity)

m_rt = art(RT ~ target * validity + Error(subject/(validity*target)), data=data) 
rt_anova = anova(m_rt)
rt_anova$peta2 = with(rt_anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))

print(rt_anova, verbose=TRUE)


### ART from aligned rank table, with aov (same results as ARTool library)
art_data <- read.csv("TemporalAttention_E5_T3_AllData_rt.art.csv")

art_data$subject <- factor(art_data$subject)
art_data$target <- factor(art_data$target)
art_data$validity <- factor(art_data$validity)

art_aov_validity <- aov(formula = ART.RT..for.validity ~ validity*target + Error(subject/(validity*target)), data = art_data)
art_aov_target <- aov(formula = ART.RT..for.target ~ validity*target + Error(subject/(validity*target)), data = art_data)
art_aov_validity.target <- aov(formula = ART.RT..for.target.validity ~ validity*target + Error(subject/(validity*target)), data = art_data)

summary(art_aov_validity)
summary(art_aov_target)
summary(art_aov_validity.target)

print(art_aov, verbose=TRUE)
