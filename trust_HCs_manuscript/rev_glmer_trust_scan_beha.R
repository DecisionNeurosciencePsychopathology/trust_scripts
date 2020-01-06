library(lme4)
library(lsmeans)
library(data.table)
library(xtable)
library(nlme)
library(car)
library(ggplot2)
library(multcomp)
library(multcompView)
library(plyr)
library(compareGroups)
library(nlme)
library(numDeriv)
library(mixed)
library(gridExtra)
library(tidyverse) #for all ggplot2 functions and subfunctions


# load data
library(readr)

#####data and demos location##########
scan_behavior = read_delim("~/Box Sync/Project Trust Game/data/processed/scan_behavior.txt",
                           "\t", escape_double = FALSE, trim_ws = TRUE)
beha_behavior = read_delim("~/Box Sync/Project Trust Game/analyses/manuscript/beha_behavior_Mar08.txt",
                          "\t", escape_double = FALSE, trim_ws = TRUE)
hall_behavior = read_delim("~/Box Sync/Project Trust Game/analyses/manuscript/hallquist_behavior.txt",
                           "\t", escape_double = FALSE, trim_ws = TRUE)

# beha_behavior = read_delim("~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/beha_behavior10302017.txt",
#                            "\t", escape_double = FALSE, trim_ws = TRUE)
# scan_behavior = read_delim("~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/scan_behavior10302017.txt",
#                            "\t", escape_double = FALSE, trim_ws = TRUE)

#load demographics
demos_scan <-read.table(file="~/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/trust_scan_demos_all.csv", sep = ',', header = TRUE)
demos_beha <-read.table(file="~/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/trust_beha_demos_all.csv", sep = ',', header = TRUE)
demos_hall <-read.table(file="~/OneDrive/Documents/Project Trust Game/analyses/manuscript hc/trust_hall_demos_all.csv", sep = ',', header = TRUE)
#demos_test <-read.table(file="~/OneDrive/Documents/Project Trust Game/analyses/manuscript hc/demos01252017.csv", sep = ',', header = TRUE)

#load coarse behavior
#coarse_beha_scan <-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/coarse_behavior_scan.csv", sep = ',', header = TRUE) 
#coarse_beha_beha <-read.table(file="~/Documents/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/coarse_behavior_beha.csv", sep = ',', header = TRUE) 

#combining beha + scan behaviors
beha_behavior$study_type = "beha"
scan_behavior$study_type = "scan"
hall_behavior$study_type = "hall"

demos_scan$study_type <- "scan"
demos_beha$study_type <- "beha"
coarse_beha_beha$study_type <- "beha"
coarse_beha_scan$study_type <- "scan"

all_behavior <- rbind(beha_behavior, scan_behavior)
all_demos <- rbind(demos_scan, demos_beha)
all_coarse <- rbind(coarse_beha_beha, coarse_beha_scan)

#looking at behavioral responses
hist(all_coarse$ratio)
ratio80.20 <- (all_coarse$ratio > 80/20 | all_coarse$ratio < 20/80)
ratio90.10 <- (all_coarse$ratio > 90/10 | all_coarse$ratio < 10/90)
count(ratio80.20) #n = 22
count(ratio90.10) #n = 21
all_coarse$ratio80.20 <- all_coarse$ratio > 80/20 | all_coarse$ratio < 20/80
all_coarse$ratio90.10 <- all_coarse$ratio > 90/10 | all_coarse$ratio < 10/90

xtabs(all_coarse$ratio80.20 ~ all_coarse$study_type) #scan n = 18
xtabs(all_coarse$ratio90.10 ~ all_coarse$study_type) #scan n = 17

#combining everything into a single table
all_coarse$subject <- all_coarse$ID
all_demos$subject <- all_demos$ID
all_data <- merge(all_behavior, all_coarse, by = "subject")
all_data <- merge(all_data, all_demos, by = "subject")

#a = data.table(beha_behavior)
#b = data.table(beha_behavior)
#b = data.table(scan_behavior)
#b = data.table(hall_behavior)
#b = data.table(hall_behavior)
#c = data.table(scan_behavior)
#b = rbind(a,b,c)

attach(b)
View(b)

#combining age w/ behavior
demos_beha$subject <- demos_beha$ID
a <- left_join(a, demos_beha, by = "subject")
any(is.na(a$AGETODAY))
demos_hall$subject <- demos_hall$ID
b <- left_join(b, demos_hall, by = "subject")
any(is.na(b$Age)) #true
unique(b$subject[is.na(b$Age)])

#lag variable: previous trustee decision
b[, pt_decision:=c(NA, t_decision[-.N]), by=trustee]

# add previous subject decision
b[, s_decision_lag1:=c(NA, s_decision[-.N]), by=trustee]
b[, s_decision_lag2:=c(NA, s_decision_lag1[-.N]), by=trustee]
b[, s_decision_lag3:=c(NA, s_decision_lag2[-.N]), by=trustee]
b[, s_decision_lag4:=c(NA, s_decision_lag3[-.N]), by=trustee]

#getting rid of uncoded trustees
b = b[b$trustee > 0]

#getting rid of missed responses
b = b[b$s_decision != 0]
#re-coding missed reponses
#b$missed = 0
#b$missed[b$s_decision == 0] = 1


#recoding variables
b$t_decision[b$t_decision==-1]=0
b$s_decision[b$s_decision==-1]=0
b$pt_decision[b$pt_decision==-1]=0
b$s_decision_lag1[b$s_decision_lag1==-1]=0


#factors
b$trustee = as.factor(b$trustee)
b$subject = as.factor(b$subject)
b$t_decision = as.factor(b$t_decision)
b$pt_decision = as.factor(b$pt_decision)
b$missed = as.factor(b$missed)

#recalculating RT
b$RT_calc = b$feedback_Onset-b$decision_Onset

#mean-centering exchange
b$mc_exchange = b$exchange - 24.5

#calculating aggregate means
#RTmeansBYsubjectBYtrustee=aggregate(b[,b$RT_calc], list(b$subject,b$trustee),mean)
RTmeansBYsubject=aggregate(b[,b$RT_calc], list(b$subject),mean)
View(RTmeansBYsubject)

attach(b)
#running glmer: w/ varying intercept and slope for subject
m0 <- glmer(s_decision ~ trustee*mc_exchange*pt_decision + (1|subject), binomial(link = "logit"), data = b,na.action = na.omit)
#m0.rev <- glmer(s_decision ~ trustee*mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit) #fails to converge
anova(m0)
summary(m0)
#for pt_decision and subject and varying slope effect of trustee within subject
m1 <- glmer(s_decision ~ trustee*mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
m2 <- glmer(s_decision ~ trustee*mc_exchange+trustee*pt_decision+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m1,m2)

#m4=best model for beha data#
m3.beh <- glmer(s_decision ~ trustee*pt_decision+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
m4.beh <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
Anova(m4.beh)
anova(m3.beh, m4.beh)
m5.beh <- glmer(s_decision ~ trustee*pt_decision+mc_exchange*pt_decision + (1|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m5.beh, m3.beh)

##the best model for scan data##
m3 <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m2,m3)
anova(m3)
Anova(m3)
summary(m3)

m4 <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m3,m4)

#best model for hallquist
m3.hall <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
m3.hall.1 <- glmer(s_decision ~ mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b[b$trustee==1,],na.action = na.omit)
m3.hall <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit,contrasts = list(trustee = "contr.sum"))
m3.hall.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b,na.action = na.omit)
ss <- getME(m3.hall.rev,c("theta","fixef"))
m3.hall.rev <- update(m3.hall.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))

contrasts(b$trustee) <- contr.treatment(levels(b$trustee),
                                      base = which(levels(b$trustee) == '2'))
rg2exch <- lsmeans(m3.hall,"mc_exchange", by = "trustee", at = list(mc_exchange = c(-24.5,  0,  24.5)))
pairs(rg2exch)

#Revisions for CABN

#reinforcement schedules
b$reward_schedule = as.factor(b$reward_schedule)
#ordering trustees
b <- arrange(b, subject, trialnum)
b$trustee_order <- 1
b$trustee_order[b$trialnum>48 & b$trialnum<=96] <- 2
b$trustee_order[b$trialnum>96 & b$trialnum<=144] <- 3
b$trustee_order[b$trialnum>144] <- 4
b$RS_order <- 1
b$RS_order[b$exchange>16 & b$exchange<=32] <- 2
b$RS_order[b$exchange>32 & b$exchange<=48] <- 3

#best model for Study 1 (or beha) is m1.rev
m0.rev <- glmer(s_decision ~ trustee*I(mc_exchange/100)*pt_decision*reward_schedule +(I(mc_exchange/100)|subject:trustee_order:RS_order), binomial(link = "logit"), data = b,nAGQ = 0,na.action = na.omit)
m1.rev <- glmer(s_decision ~ trustee*pt_decision*reward_schedule +pt_decision*I(mc_exchange/100) +(I(mc_exchange/100)|subject:trustee_order:RS_order), binomial(link = "logit"), data = b,nAGQ = 0,na.action = na.omit)
m1.rev <- glmer(s_decision ~ trustee*pt_decision*reward_schedule +pt_decision*I(mc_exchange/100) +(I(mc_exchange/100)|subject:trustee_order:RS_order), binomial(link = "logit"), data = b,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),na.action = na.omit)
#m1a.rev <- glmer(s_decision ~ trustee*pt_decision*reward_schedule +pt_decision*I(mc_exchange/100) +(1+pt_decision|subject:trustee_order:RS_order), binomial(link = "logit"), data = b,nAGQ = 0,na.action = na.omit)
m2.rev <- glmer(s_decision ~ trustee*pt_decision+trustee*reward_schedule+pt_decision*reward_schedule +pt_decision*I(mc_exchange/100) +(I(mc_exchange/100)|subject:trustee_order:RS_order), binomial(link = "logit"), data = b,nAGQ = 0,na.action = na.omit)

#better model w/out mc_exchange which is conflated with RS order: best model for Study 1 (beha) = m0.rev
m0.rev <- glmer(s_decision ~ trustee*pt_decision*reward_schedule +(1|subject:trustee_order:RS_order), binomial(link = "logit"), data = b,nAGQ = 0,na.action = na.omit)
m0.rev <- glmer(s_decision ~ trustee*pt_decision*reward_schedule +(1|subject:trustee_order:RS_order), binomial(link = "logit"), data = b,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),na.action = na.omit)
m1.rev <- glmer(s_decision ~ trustee*pt_decision+trustee*reward_schedule+ reward_schedule*pt_decision+(1|subject:trustee_order:RS_order), binomial(link = "logit"), data = b,nAGQ = 0,na.action = na.omit)

#comparison of best model from before with RS added; significant difference betwen AIC, indicating that the model with reward_schedule explains more variance (Study 1: best model m6b.rev)
#Study 1 = beha
m3.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
m4.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + reward_schedule+(1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m4.rev,c("theta","fixef"))
m4a.rev <- update(m4.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
m5.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + reward_schedule*pt_decision+(1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m5.rev,c("theta","fixef"))
m5a.rev <- update(m5.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
m6.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + reward_schedule*pt_decision*trustee+(1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m6.rev,c("theta","fixef"))
m6a.rev <- update(m6.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
ss <- getME(m6a.rev,c("theta","fixef"))
m6b.rev <- update(m6a.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))

#Looking at s_decision_lag1 or PSD == 1 (INVESTED), R2, Q12
m14.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 1,],na.action = na.omit)
ss <- getME(m14.rev,c("theta","fixef"))
m14.rev <- update(m14.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
summary(m14.rev)
#note that the pt_decision is not significantly different when random effects are defined as below (I don't know what to make of this)
m14a.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b[b$s_decision_lag1 == 1,],na.action = na.omit)
m14b.rev <- glmer(s_decision ~ trustee+mc_exchange+pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b[b$s_decision_lag1 == 1,],na.action = na.omit)
anova(m14a.rev, m14b.rev) #ns

#R2, Q14: PSD == 0 (KEPT)
m17.rev <- glmer(s_decision ~ trustee+pt_decision + pt_decision*reward_schedule+(1+pt_decision|subject:trustee_order:RS_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
m18.rev <- glmer(s_decision ~ trustee+pt_decision + pt_decision+(1+pt_decision|subject:trustee_order:RS_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
m19.rev <- glmer(s_decision ~ trustee+pt_decision*mc_exchange + pt_decision*reward_schedule+(1+pt_decision|subject), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m19.rev,c("theta","fixef"))
m19.rev <- update(m19.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4))) #same result


#Study 2 = hallquist
m9.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m9.rev,c("theta","fixef"))
m9a.rev <- update(m9.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
m10.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + reward_schedule +(1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m10.rev,c("theta","fixef"))
m10a.rev <- update(m10.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))

#Looking at s_decision_lag1 or PSD == 1 (INVESTED), R2, Q12
m14.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 1,],na.action = na.omit)
ss <- getME(m14.rev,c("theta","fixef"))
m14.rev <- update(m14.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
summary(m14.rev)
m14a.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b[b$s_decision_lag1 == 1,],na.action = na.omit)
summary(m14a.rev) #same result for pt_decision1 as with the other model

#Looking at s_decision_lag1 or PSD == 0 (KEPT), R2, Q13
m15.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0,],na.action = na.omit)
ss <- getME(m15.rev,c("theta","fixef"))
m15.rev <- update(m15.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
summary(m15.rev)
m16.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0,],na.action = na.omit)
summary(m16.rev)

#R2, Q14: PSD == 0 (KEPT)
m17.rev <- glmer(s_decision ~ trustee+pt_decision + pt_decision*reward_schedule+(1+pt_decision|subject:trustee_order:RS_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m17.rev,c("theta","fixef"))
m17.rev <- update(m17.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
m18.rev <- glmer(s_decision ~ trustee+pt_decision + pt_decision+(1+pt_decision|subject:trustee_order:RS_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
m19.rev <- glmer(s_decision ~ trustee+pt_decision*mc_exchange + pt_decision*reward_schedule+(1+pt_decision|subject), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m19.rev,c("theta","fixef"))
m19.rev <- update(m19.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4))) #same result

#Study 3 = scan
m7.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m7.rev,c("theta","fixef"))
m7a.rev <- update(m7.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
m8.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + reward_schedule+(1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m8.rev,c("theta","fixef"))
m8a.rev <- update(m8.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))

#test <- glmer(s_decision ~ reward_schedule+(1|subject), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)

#Looking at s_decision_lag1 or PSD == 1 (INVESTED), R2, Q12
m14.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 1,],na.action = na.omit)
ss <- getME(m14.rev,c("theta","fixef"))
m14.rev <- update(m14.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
summary(m14.rev)
m14a.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b[b$s_decision_lag1 == 1,],na.action = na.omit)
#same result for pt_decision as for m14.rev

#R2, Q14: PSD == 0 (KEPT)
m17.rev <- glmer(s_decision ~ trustee+pt_decision + pt_decision*reward_schedule+(1+pt_decision|subject:trustee_order:RS_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m17.rev,c("theta","fixef"))
m17.rev <- update(m17.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
m18.rev <- glmer(s_decision ~ trustee+pt_decision + pt_decision+(1+pt_decision|subject:trustee_order:RS_order), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
m19.rev <- glmer(s_decision ~ trustee*mc_exchange+pt_decision*mc_exchange + pt_decision*reward_schedule+(1+pt_decision|subject), binomial(link = "logit"), data = b[b$s_decision_lag1 == 0 & b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m19.rev,c("theta","fixef"))
m19.rev <- update(m19.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4))) #same result

#scan + hallquist combined: m12a.rev best for combined set
m11.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m11.rev,c("theta","fixef"))
m11a.rev <- update(m11.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
m12.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m12.rev,c("theta","fixef"))
m12a.rev <- update(m12.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4))) 
m13.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + reward_schedule+(1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b[b$reward_schedule!=50,],na.action = na.omit)
ss <- getME(m13.rev,c("theta","fixef"))
m13a.rev <- update(m13.rev, start = ss, control=glmerControl(optCtrl=list(maxfun=2e4)))

#R2, Q15: Studies 1 and 3 combined
m20.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
m21.rev <- glmer(s_decision ~ trustee+mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m20.rev, m21.rev) #p =0.02939
Anova(m20.rev,"3")
#m22.rev <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + (1+pt_decision|subject:trustee_order), binomial(link = "logit"), data = b,na.action = na.omit)
#trustee is not significant in this model, likely because trustee_order is conflated with trustee effects. I don't think that the glmer can effectively differentiate between the nested
#trustee_order slope within subject and the trustee effects.

#GRAPHS
#trustee
pd = position_dodge(.05)
lsm<- emmeans::lsmeans(m20.rev, ~trustee)
pairs(lsm)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey", interaction = "TRUE")
em$trustee = revalue(em$trustee , c("1" = "good","2"="bad","3"="neutral","4" = "computer"))
ggplot(data = em, aes(x=trustee, y=lsmean, color = trustee)) +
  geom_errorbar(aes(ymin = lsmean-SE, ymax = lsmean+SE), width = .2) +
  geom_point(position=pd)+
  xlab("Trustee Type") + ylab("Estimated Likelihood of Participant's Investment")

pd = position_dodge(.05)
lsm<- emmeans::lsmeans(m4.beh, ~trustee)
pairs(lsm)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey", interaction = "TRUE")
em$trustee = revalue(em$trustee , c("1" = "good","2"="bad","3"="neutral","4" = "computer"))
ggplot(data = em, aes(x=trustee, y=lsmean, color = trustee)) +
  geom_errorbar(aes(ymin = lsmean-SE, ymax = lsmean+SE), width = .2) +
  geom_point(position=pd)+
  xlab("Trustee Type") + ylab("Estimated Likelihood of Participant's Investment")


#trustee x pt_decision
pd = position_dodge(.05)
lsm<- emmeans::lsmeans(m0.rev,~pt_decision|trustee)
contrast(lsm, method = "pairwise", adjust = "tukey",interaction = TRUE)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey", interaction = "TRUE")
em$trustee = revalue(em$trustee , c("1" = "good","2"="bad","3"="neutral","4" = "computer"))
em$pt_decision = revalue(em$pt_decision , c("0" = "KEPT","1"="SHARED"))
ggplot(data = em, aes(x=pt_decision, y=lsmean, color = trustee)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .2) +
  geom_point(position=pd)+
  facet_wrap(~trustee, ncol=4)+
  xlab("Trustee's Decision on Previous Trial") + ylab("Estimated Likelihood of Participant's Investment")

#trustee x reward schedule
pd = position_dodge(.05)
lsm<- emmeans::lsmeans(m6b.rev,~reward_schedule|trustee)
contrast(lsm, method = "pairwise", adjust = "tukey",interaction = TRUE)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey", interaction = "TRUE")
em$trustee = revalue(em$trustee , c("1" = "good","2"="bad","3"="neutral","4" = "computer"))
ggplot(data = em, aes(x=reward_schedule, y=lsmean, color = trustee)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .2) +
  geom_point(position=pd)+
  facet_wrap(~trustee, ncol=4)+
  xlab("Reward Schedule") + ylab("Estimated Likelihood of Participant's Investment")

# trustee*mc_exchange
lsm <- emmeans::lsmeans(m3, "mc_exchange", by = "trustee",at = list(mc_exchange = c(-24.5,  0,  24.5)))
lsm <- emmeans::lsmeans(m3.hall, "mc_exchange", by = "trustee",at = list(mc_exchange = c(-24.5,  0,  24.5)))
contrast(lsm, method = "pairwise", adjust = "tukey",interaction = TRUE)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey", interaction = "TRUE")
em$trustee = revalue(em$trustee , c("1" = "good","2"="bad","3"="neutral","4" = "computer"))
ggplot(data = em, aes(x=mc_exchange, y=lsmean, color = trustee)) +
  geom_errorbar(aes(ymin = lsmean-SE, ymax = lsmean+SE), width = .2) +
  geom_point(position=pd)+
  facet_wrap(~trustee, ncol=4)+
  xlab("Number of Exchanges with Trustee") + ylab("Estimated Likelihood of Participant's Investment")

#exchanges x PTD
pd = position_dodge(.05)
lsm<- emmeans::lsmeans(m6b.rev,~mc_exchange|pt_decision,at = list(mc_exchange = c(-95,  0,  95)))
contrast(lsm, method = "pairwise", adjust = "tukey",interaction = TRUE)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey")
em$trustee = revalue(em$trustee , c("1" = "good","2"="bad","3"="neutral","4" = "computer"))
ggplot(data = em, aes(x=mc_trialnum, y=lsmean, color = trustee)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position=pd, width = .2) +
  geom_point(position=pd)+
  facet_wrap(~trustee, ncol=4)+
  scale_y_reverse()+
  xlab("Practice (Mean-Centered Trial Number)") + ylab("Estimated RT (1000/x)")

#reward schedule x pt_decision
pd = position_dodge(.05)
lsm<- emmeans::lsmeans(m0.rev,~pt_decision|reward_schedule)
contrast(lsm, method = "pairwise", adjust = "tukey",interaction = TRUE)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey", interaction = "TRUE")
em$pt_decision = revalue(em$pt_decision , c("0" = "KEPT","1"="SHARED"))
ggplot(data = em, aes(x=pt_decision, y=lsmean, color = reward_schedule)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .2) +
  geom_point(position=pd)+
  facet_wrap(~reward_schedule, ncol=4)+
  xlab("Trustee's Decision on Previous Trial") + ylab("Estimated Likelihood of Participant's Investment")

#trustee x reward schedule x pt_decision
lsm<- emmeans::lsmeans(m0.rev,~pt_decision|reward_schedule*trustee)
lsm<- emmeans::lsmeans(m6b.rev,~pt_decision|reward_schedule*trustee)
contrast(lsm, method = "pairwise", adjust = "tukey",interaction = TRUE)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey", interaction = "TRUE")
em$trustee = revalue(em$trustee , c("1" = "good","2"="bad","3"="neutral","4" = "computer"))
em$pt_decision = revalue(em$pt_decision , c("0" = "KEPT","1"="SHARED"))
ggplot(data = em, aes(x=pt_decision, y=lsmean, color = reward_schedule)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .2) +
  geom_point(position=pd)+
  facet_wrap(reward_schedule~trustee, ncol=4)+
  xlab("Trustee's Decision on Previous Trial") + ylab("Estimated Likelihood of Participant's Investment")


#establishing a reference grid
likelihood.m3 <- ref.grid(m3)

#marginal means of trustree, plus contrast
trustee.lsm <- lsmeans(likelihood.m3, "trustee")
trustee.sum <- summary(trustee.lsm, infer = c(TRUE,TRUE), level = 0.95, adjust = "bon")
contrast(trustee.lsm, method = "pairwise", adjust ="bon")

# plot main effect of trustee
trustee.plot1 <- plot(lsmeans(likelihood.m3, "trustee"))
trustee.plot2 <- plot(trustee.lsm, type ~ trustee, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Trustee Type")
#axis(side = "bottom",labels = c("good","bad","neutral","computer"), at=1:4)

#rg2trust <- lsmeans(m3,by = "trustee")
#pairs(rg2trust)
#plot(rg2trust, type ~ trustee, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Trustee Type")

# plot trustee*mc_exchange
trusteeXexchange.lsm <- lsmeans(likelihood.m3, "mc_exchange", by = "trustee")
trusteeXexchange.sum <- summary(trusteeXexchange.lsm, infer = c(TRUE,TRUE), level = 0.95, adjust = "bon")

rg2exch <- lsmeans(m3,"mc_exchange", by = "trustee", at = list(mc_exchange = c(-24.5,  0,  24.5)))
plot(rg2exch, type ~ trustee, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Number of Exchanges")
pairs(rg2exch)
rg2exchB <- lsmeans(m3,"mc_exchange", by = "trustee", at = list(mc_exchange = c(-24.5, 24.5)))
pairs(rg2exchB)

# plot pt_decision*mc_exchange
rg3exch <- lsmeans(m3,"pt_decision", by = "mc_exchange", at = list(mc_exchange = c(-24.5,  0,  24.5)))
pt_decisionXmc_exchange.plot <- plot(rg3exch, type ~ pt_decision, horiz=F, ylab = "Probability of sharing (predicted)", xlab = "Previous trustee decision")
pairs(rg3exch)
pt_decisionXmc_exchange.plot$lwd=2
pt_decisionXmc_exchange.plot$lwd

# testing IF p(s_decision = SHARE) on trials following trials w/ s_decision = KEEP 
# increased FOR t_decision = SHARE vs KEEP on the current trial
#g1 <- glmer(ns_decision ~ trustee*mc_exchange*pt_decision + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)


## test pt_decision*ps_decision
# first, test model with lags

m3lag <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + pt_decision*s_decision_lag1 + s_decision_lag2 + s_decision_lag3 + (1+pt_decision|subject), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m3lag)
summary(m3lag)

anova(m3,m3lag)

## ad: this looks currently (1/19/17) like the best-fitting model
# lag means that lagged choice is included and t means that trustee is nested within subject on the random side
m3lagt <- glmer(s_decision ~ trustee*mc_exchange+mc_exchange*pt_decision + pt_decision*s_decision_lag1 + (1+pt_decision|subject) + (1|subject:trustee), binomial(link = "logit"), data = b,na.action = na.omit)
anova(m3lagt)
summary(m3lagt)

##analyzing RTs##
b$RT_calc_ln = log(b$RT_calc)
boxplot(split(b$RT_calc_ln,b$subject))

attach(b)
m0.RT <- lme(RT_calc_ln ~ trustee*mc_exchange*pt_decision*s_decision, random = (1|subject), data = b,na.action = na.omit)
m1.RT <- lme(RT_calc_ln ~ trustee*pt_decision*mc_exchange*s_decision + (1+pt_decision|subject), data = b,na.action = na.omit)
m2.RT <- lme(RT_calc_ln ~ trustee*pt_decision*mc_exchange*s_decision + (1+pt_decision|subject)+(1+s_decision|subject), data = b,na.action = na.omit)
m3.RT <- lme(RT_calc_ln ~ trustee*pt_decision*mc_exchange*s_decision + (1+s_decision|subject), data = b,na.action = na.omit)

##analyzing missed responses###
missedBYsubject=aggregate(b[,b$missed], list(b$subject),sum)
m0.missed <-glmer(missed ~ trustee*mc_exchange + (1|subject), binomial(link = "logit"), data = b,na.action = na.omit)
