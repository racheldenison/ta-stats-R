# temporal_attention_E3_adjust_nonparametric_anova.R

setwd("/Users/rachel/Documents/NYU/Projects/Temporal_Attention/XLS")
data <- read.csv("TemporalAttention_E3_adjust_AllData_swapNoBiasMAP.csv")

# http://depts.washington.edu/aimgroup/proj/art/
library(ARTool)

data$subject <- factor(data$subject)
data$T1T2 <- factor(data$T1T2)
data$validity <- factor(data$validity)

m_sd = art(sd ~ T1T2 * validity + Error(subject/(validity*T1T2)), data=data) 
sd_anova = anova(m_sd)
sd_anova$peta2 = with(sd_anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))

m_g = art(g ~ T1T2 * validity + Error(subject/(validity*T1T2)), data=data) 
g_anova = anova(m_g)
g_anova$peta2 = with(g_anova, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))

print(sd_anova, verbose=TRUE)
print(g_anova, verbose=TRUE)


### ART from table, with aov
art_data <- read.csv("TemporalAttention_E3_adjust_AllData_mixtureNoBiasMAP_sd.art.csv")

art_data$subject <- factor(art_data$subject)
art_data$T1T2 <- factor(art_data$T1T2)
art_data$validity <- factor(art_data$validity)

art_aov_validity <- aov(formula = ART.sd..for.validity ~ validity*T1T2 + Error(subject/(validity*T1T2)), data = art_data)
art_aov_T1T2 <- aov(formula = ART.sd..for.T1T2 ~ validity*T1T2 + Error(subject/(validity*T1T2)), data = art_data)
art_aov_validity.T1T2 <- aov(formula = ART.sd..for.T1T2.validity ~ validity*T1T2 + Error(subject/(validity*T1T2)), data = art_data)

summary(art_aov_validity)
summary(art_aov_T1T2)
summary(art_aov_validity.T1T2)

art_aov_validity2 <- aov(formula = ART.sd..for.validity ~ validity*T1T2 + Error(subject/(validity*T1T2)), data = art_data)


### ART example
higgins_data <-read.csv("Higgins1990-Table5.art.csv")

higgins_data$Tray <- factor(higgins_data$Tray)
higgins_data$Fertilizer <- factor(higgins_data$Fertilizer)
higgins_data$Moisture <- factor(higgins_data$Moisture)

m = art(DryMatter ~ Fertilizer * Moisture + Error(Tray/Fertilizer) + Moisture, data=higgins_data)
anova(m)
#m = art(DryMatter ~ Fertilizer * Moisture + (1|Tray), data=higgins_data)

aov_m <- aov(formula = ART.DryMatter..for.Fertilizer ~ Fertilizer*Moisture + Error(Tray/Fertilizer) + Moisture, data = higgins_data)
summary(aov_m)
aov_m <- aov(formula = ART.DryMatter..for.Moisture ~ Fertilizer*Moisture + Error(Tray), data = higgins_data)
aov_m <- aov(formula = ART.DryMatter..for.Moisture.Fertilizer ~ Fertilizer*Moisture + Error(Tray), data = higgins_data)


### T1 and T2 separately
data_T1 <- subset(data, T1T2==1)

m_sd_T1 = art(sd ~ validity + (1|subject), data=data_T1)
sd_anova_T1 = anova(m_sd_T1)

m_g_T1 = art(g ~ validity + (1|subject), data=data_T1)
g_anova_T1 = anova(m_g_T1)

data_T2 <- subset(data, T1T2==2)

m_sd_T2 = art(sd ~ validity + (1|subject), data=data_T2)
sd_anova_T2 = anova(m_sd_T2)

m_g_T2 = art(g ~ validity + (1|subject), data=data_T2)
g_anova_T2 = anova(m_g_T2)


### pairwise comparisons, p-value adjustment
t1sd <- c(0.02686, 0.00928, 0.85010)
p.adjust(t1sd, method = 'holm')

t2sd <- c(0.17627, 0.56934, 0.12939)
p.adjust(t2sd, method = 'holm')

t1g <- c(0.67725, 0.67725, 0.96973)
p.adjust(t1g, method = 'holm')

t2g <- c(0.20361, 0.67725, 0.79102)
p.adjust(t2g, method = 'holm')

# circular sd (raw distributions)
t2circsd <- c(0.00472, 0.93244, 0.00932)
p.adjust(t2circsd, method = 'holm')

# session 1
t2g_sess1 <- c(0.02686, 0.06396, 0.09229)
p.adjust(t2g_sess1, method = 'holm')

t1sd_sess1 <- c(0.00488, 0.10986, 0.26611)
p.adjust(t1sd_sess1, method = 'holm')
