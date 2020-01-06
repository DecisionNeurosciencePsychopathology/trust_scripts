library(lme4)
library(lsmeans)
library(data.table)
library(nlme)
library(xtable)
library(ggplot2)
library(multcomp)
library(multcompView)
library(plyr)
# load data
library(readr)

#Merge all the trust betas into one dataset

#Load in the data (betas from PE maps produced by mean-centered but not z-scored PE regressors)
df_null<-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/betas/new_betas_null.txt", sep = '\t', header = TRUE)
df_regret<-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/betas/new_betas_regret.txt", sep = '\t', header = TRUE)
df_trustee<-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/betas/new_betas_trustee.txt", sep = '\t', header = TRUE)
df_policy<-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/betas/new_betas_policy.txt", sep = '\t', header = TRUE)

#betas from z-scores
df_null<-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/betas/betas_zscored_real_Pes.txt", sep = '\t', header = TRUE)
df_regret<-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/betas/betas_zscored_counter_hybrid_regret_Pes.txt", sep = '\t', header = TRUE)
df_trustee<-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/betas/betas_zscored_counter_trustee_Pes.txt", sep = '\t', header = TRUE)
df_policy<-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/betas/betas_zscored_subject_counter_Pes.txt", sep = '\t', header = TRUE)

#load ages
ages <- read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript/ids_ages.txt", sep = '\t', header = TRUE)

#Set model variables
models <- c("null", "regret", "trustee","policy")

#Remove all nans
dfList <- list('df_null','df_regret', 'df_trustee','df_policy')
l<-lapply(dfList,function(x){
    x <- get(x)
    x <- na.omit(x)
    return(x)
}
)

for ( i in 1:length(dfList)) {
    l[[i]]$model <- rep(models[i],nrow(l[[i]])) #Assign model names
    assign(dfList[[i]], l[[i]])
}

df_policy <- merge(df_policy, ages, by = "ID")
df_null <- merge(df_null, ages, by = "ID")
df_regret <- merge(df_regret, ages, by = "ID")
df_trustee <- merge(df_trustee, ages, by = "ID")

#Make wide and long version of data
long_trust_betas <- rbind(df_policy,df_null,df_regret,df_trustee)
long_trust_betas <- melt(long_trust_betas, id.vars = c("ID","model","age_consent", "age_today"))
colnames(long_trust_betas)[5] <- "ROI"
colnames(long_trust_betas)[6] <- "beta"


wide_trust_betas <- dcast(long_trust_betas, ID + age_consent + age_today ~ model + ROI, value.var = "beta")

#saving raw data in long and wide formats
save(long_trust_betas, file = "long_format_beta_coef.RData")
save(wide_trust_betas, file = "wide_format_beta_coef.RData")

##saving the z-scores in long and wide format
save(long_trust_betas, file = "long_format_zbeta_coef.RData")
save(wide_trust_betas, file = "wide_format_zbeta_coef.RData")


m0.lmer <- lmer(beta ~ model + ROI + (1|ID), long_trust_betas)
m0b.lmer <- lmer(beta ~ model + ROI + (1|ID) + (1|ID:ROI), long_trust_betas) ##this one## Random effects: intercept for ID and intercept of ROI within ID (nested term) 
m0c.lmer <- lmer(beta ~ model * ROI + (1|ID) + (1|ID:ROI), long_trust_betas) 
m0d.lmer <- lmer(beta ~ model + ROI + (1+ROI|ID), long_trust_betas) 
lsmeans.model.m0b <- lsmeans(m0b.lmer, "model")
lsmeans.modelXroi.m0b <- lsmeans(m0b.lmer, "model", by = "ROI")
contrast(lsmeans.model.m0b, method = "pairwise", adjust = "tukey")

Anova(m0.lmer)
summary(m0.lmer)

lsmeans.model.m0 <- lsmeans(m0.lmer,"model")
lsmeans.ROI.m0 <- lsmeans(m0.lmer,"ROI")
lsmeans.modelBYroi.m0 <- lsmeans(m0.lmer, "model", by = "ROI")


plot_actual_betas <- ggplot(aes(x=ROI, y=beta, color = model),data = long_trust_betas) +
    geom_boxplot()

##re-formatting lsmeans object for plotting lsmeans using ggplot, cld uses library(plyr)
em_betas.m0b = cld(lsmeans.modelBYroi.m0b, alpha = 0.05, Letters = letters, adjust = "tukey")
em_betas.m0b$.group = gsub(" ", "", em_betas.m0b$.group)

#em_betas.m0b <- summary(lsmeans.modelXroi.m0b) #alternative way of transforming the lsmeans object to be used for ggplot
em_betas.m0b$model = revalue(em_betas.m0b$model, c("null" = "actual"))
em_betas.m0b$ROI = revalue(em_betas.m0b$ROI, c("Ventral.Striatum" = "Ventral Striatum", "Left.Caudate" = "Left Caudate", "Frontal.Operculum" = "Operculum"))
pd = position_dodge(0.5)
plot_est_betas.a <- ggplot(data=em_betas.m0b, aes(x = ROI, y=lsmean, color = model)) +
                    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), position = pd, width = .2) +
                    geom_point(position = pd)
em_betas.m0b = cld(lsmeans.model.m0b, alpha = 0.05, Letters = letters, adjust = "tukey")
em_betas.m0b$model = revalue(em_betas.m0b$model, c("null" = "actual"))
em_betas.m0b$model = factor(em_betas.m0b$model, levels = c("actual", "regret", "trustee", "policy"))
plot_est_betas.b <- ggplot(data=em_betas.m0b, aes(x = model, y = lsmean)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),width = .1) +
    geom_point() +
    theme_bw() +
    xlab("Model Type") + ylab("Estimated Mean of Beta Coefficients") +
    theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11))
save(plot_est_betas.b, file = "plot_est_betas.b")
    
##################################not reported in manuscript###########################################    
##older analyses, w/ and w/out outliers, etc.##
contrast(lsmeans.model.m0, method = "pairwise", adjust ="bon")
contrast(lsmeans.ROI.m0, method = "pairwise", adjust ="bon")
contrast(lsmeans.modelBYroi.m0, method = "eff", adjust = "tukey",interaction = FALSE)

plot(lsmeans.model.m0)
plot(lsmeans.ROI.m0)

m1.lmer <- lmer(beta ~ model*ROI + (1|ID), long_trust_betas)
anova(m0.lmer, m1.lmer) #ns, p = 0.9214

m2.lmer <- lmer(beta ~ model + ROI + age_today + (1|ID), long_trust_betas)
lsmeans.model.m2 <- lsmeans(m2.lmer,"model")
plot(lsmeans.model.m2)
contrast(lsmeans.model.m2, method = "pairwise", adjust ="tukey")

##wout occipital
data_ltb = subset(long_trust_betas, ROI != "Occipital")
m3a.lmer <- lmer(beta ~ model + ROI + age_today + (1|ID), data_ltb)
Anova(m3a.lmer)
lsmeans.model.m3a <- lsmeans(m3a.lmer,"model")
plot(lsmeans.model.m3a)
#contrast(lsmeans.model.m3a, method = "pairwise", adjust ="tukey")
#lsmeans.modelBYroi.m3a<- lsmeans(m3a.lmer, "model", by = "ROI")

##simple box plot of actual data
box_parameters <- boxplot(beta ~ model+ROI,data = long_trust_betas)
    
wout_outliers <- wide_trust_betas[,2:ncol(wide_trust_betas)]
for (i in 1:ncol(box_parameters$stats))
{
    wout_outliers[wout_outliers[,i] < box_parameters$stats[1] | wout_outliers[,i] > box_parameters$stats[5],i] <- NA
}
wout_outliers <- cbind(wide_trust_betas[,1],wout_outliers)
colnames(wout_outliers)[1] <- "ID"
l1 <- wout_outliers[1:6]
l1$model <- "null"
colnames(l1)[2] <- "Ventral.Striatum"
colnames(l1)[3] <- "Frontal.Operculum"
colnames(l1)[4] <- "Thalamus"
colnames(l1)[5] <- "Occipital"
colnames(l1)[6] <- "Left.Caudate"

l2 <- wout_outliers[c(1,7:11)]
l2$model <- "policy"
colnames(l2)[2] <- "Ventral.Striatum"
colnames(l2)[3] <- "Frontal.Operculum"
colnames(l2)[4] <- "Thalamus"
colnames(l2)[5] <- "Occipital"
colnames(l2)[6] <- "Left.Caudate"

l3 <- wout_outliers[c(1,12:16)]
l3$model <- "regret"
colnames(l3)[2] <- "Ventral.Striatum"
colnames(l3)[3] <- "Frontal.Operculum"
colnames(l3)[4] <- "Thalamus"
colnames(l3)[5] <- "Occipital"
colnames(l3)[6] <- "Left.Caudate"

l4 <- wout_outliers[c(1,17:21)]
l4$model <- "trustee"
colnames(l4)[2] <- "Ventral.Striatum"
colnames(l4)[3] <- "Frontal.Operculum"
colnames(l4)[4] <- "Thalamus"
colnames(l4)[5] <- "Occipital"
colnames(l4)[6] <- "Left.Caudate"

l_all <- rbind(l1, l2, l3, l4)
l_wout_outliers <- melt(l_all, id.vars = c("ID","model"))
colnames(l_wout_outliers)[3] <- "ROI"
colnames(l_wout_outliers)[4] <- "beta"

m3.lmer <- lmer(beta ~ model + ROI + (1|ID), l_wout_outliers)
m4.lmer <- lmer(beta ~ model*ROI + (1|ID), l_wout_outliers) #ns, p = 0.8834
Anova(m3.lmer)
summary(m3.lmer)

lsmeans.model.m3 <- lsmeans(m3.lmer,"model")
lsmeans.ROI.m3 <- lsmeans(m3.lmer,"ROI")

plot(lsmeans.model.m3)
plot(lsmeans.ROI.m3)

contrast(lsmeans.model.m3, method = "pairwise", adjust ="tukey")
contrast(lsmeans.ROI.m3, method = "pairwise", adjust ="tukey")
#############################################################################    
