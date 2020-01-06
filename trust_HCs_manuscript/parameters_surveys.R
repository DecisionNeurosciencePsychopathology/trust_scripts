library(lme4)
library(lsmeans)
library(data.table)
library(nlme)
library(car)
library(stats4)
library(readr)
library(reshape2)
library(languageR)
library(lattice)
library(plyr)
library(ggplot2)
library(multcomp)
library(multcompView)

#reading in data
data = read_delim("~/Box Sync/Project Trust Game/analyses/manuscript/surveys and parameters uncensored.txt",
           "\t", escape_double = FALSE, trim_ws = TRUE)
data$kc[data$kc == 0] = NA
data$deltaKGKB = data$kg-data$kb
data$deltaRGRB = data$rTGpR-data$rTtBpR

#boxplots of parameters by condition
box_parameters <- boxplot(data$beta, data$ks, data$kg, data$kb, data$kn, data$theta, col = "deepskyblue", names = c("beta","ks","kg","kb","kn","theta"))

#histograms
data.ks.hist <- hist(ks)
data.kg.hist <- hist(kg) 
data.kb.hist <- hist(kb)
data.kn.hist <- hist(kn)
data.beta.hist <- hist(beta)
data.theta.hist <- hist(theta)

#scatterplots
parXpar.scatter = scatterplot(deltaKGKB,deltaRGRB);
parXpar <- pairs(~ks+kg+kb+kn, data = data)
ksXrating <- pairs(~ks+rTGpR+rTtBpR+rTtNpR, data = data);

#factors
data$ID = as.factor(data$ID)

#reshaping the dataset
dfparams <- melt(data[,c(1:6,8)],id.vars = c("ID","beta","ks","theta"),value.name = "kappa")
names(dfparams)[names(dfparams)=="variable"] <- "trustee"
dfparams <- na.omit(dfparams)
dfparams$trustee <- revalue(dfparams$trustee, c("kg" ="good","kb"="bad","kn"="neutral"))

dfratings <- melt(data[,c(1,10:12)],id.vars = c("ID"),value.name = "rating")
names(dfratings)[names(dfratings)=="variable"] <- "trustee"
dfratings$trustee <- revalue(dfratings$trustee, c("rTGpR" ="good","rTtBpR"="bad","rTtNpR"="neutral"))

df <- merge(dfparams,dfratings)

#checking factors
df$ID = as.factor(df$ID)
df$trustee = as.factor(df$trustee)

#Predicting rating with parameters + trustee type  
m0.lmer = lmer(rating ~ trustee+ks+kappa +(1|ID), df)
m1.lmer = lmer(rating ~ trustee*ks+kappa+ (1|ID), df)
m2.lmer = lmer(rating ~ trustee+kappa+ (1|ID), df) #this one (not used in manuscript)
m3.lmer = lmer(rating ~ trustee*kappa+ (1|ID), df)
anova(m2.lmer, m3.lmer) #ns

#Predicting condition-level parameter with trustee type, subject-level parameter and rating
m4.lmer = lmer(kappa ~ trustee + ks + rating + (1|ID), df) 
m5.lmer = lmer(kappa ~ trustee + ks  + (1|ID), df) 
m6.lmer = lmer(kappa ~ rating + ks  + (1|ID), df)#this one reported in manuscript

ls4 <- lsmeans(m4.lmer,"trustee")
plot(ls4)

ls6 <- lsmeans(m6.lmer,"rating", at = list(rating = c(1,4,7))) #this one reported in manuscript
ls7 <- lsmeans(m4.lmer,"rating", at = list(rating = c(1,4,7)))
plot(ls6)
plot(ls7)

em.plot <- plot(ls6, xlab = "Estimated Means of Trustee-level Bias Parameter", ylab = "Low, Mid-, and High Trustworthiness Ratings", col = "black")
#lsmip(m3.lmer, trustee ~ kappa)

likelihood.m2 <- ref.grid(m2.lmer)
trustee.lsm <- lsmeans(likelihood.m2, "trustee")
trustee.sum <- summary(trustee.lsm, infer = c(TRUE,TRUE), level = 0.95, adjust = "bon")
contrast(trustee.lsm, method = "pairwise", adjust ="bon")

Anova(m2.lmer)

#Checking that the parameters are significantly different from each other (also reported in manuscript)
m7.lmer <- lmer(kappa ~ trustee + (1|ID),df)
kappa.lsmeans <- lsmeans(m7.lmer, "trustee")
contrast(kappa.lsmeans, method = "pairwise", adjust ="tukey")

#better plotting of lsmeans
em_bias = cld(ls6, alpha = 0.05, Letters = letters, adjust = "tukey")
em_bias$rating <- factor(em_bias$rating)
plot_est_bias <- ggplot(data=em_bias, aes(x = rating, y = lsmean)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),width = .1) +
    geom_point() +
    theme_bw() +
    xlab("Trustworthiness Rating") + ylab("Estimated Mean of Trustee-Level Bias") +
    scale_x_discrete(breaks=c("1","4","7"),labels = c("Low","Mid", "High")) +
    theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) 
    