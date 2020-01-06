library(readr) #read data
library(ggplot2)
library(gridExtra)
library(tidyverse) #for all ggplot2 functions and subfunctions
library(dplyr)
library(plyr)
library(compareGroups)
library(lme4)
library(lmerTest)
#library(lsmeans)
library(emmeans)
library(data.table)
library(xtable)
library(nlme)
library(car)
library(multcomp)
library(multcompView)
library(numDeriv)
library(mixed)
library(MASS)
library(effects)
library(readr)
library(VIM)
library(mice)
library(multcompView)
library(stargazer)

df_policy<-read.table(file="~/OneDrive/Documents/Project Trust Game/analyses/manuscript hc/betas/betas_zscored_subject_counter_Pes.txt", sep = '\t', header = TRUE)
betas_ROIS = read_delim("PE_ROIbetas_dec_feed.txt",
                           "\ ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
betas_contrasts_ROIs = read_delim("decXfeed_betas.txt",
                                  "\ ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
betas_ROIS2 = read_delim("PE_ROIbetas_dec_feed_wTAUSQR.txt",
                        " ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
betas_contrasts_ROIS2 = read_delim("decXfeed_betas_wTAUSQR.txt",
                                  "\ ", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

clusters <- kmeans(betas_ROIS[,1:3], 5)
clusters <- kmeans(betas_ROIS2[,1:3], 5)
betas_ROIS$cluster <- as.factor(clusters$cluster)
betas_ROIS2$cluster <- as.factor(clusters$cluster)
names(betas_ROIS2) [7] <- "PE-decision"
names(betas_ROIS2) [8] <- "dec_tau2"
names(betas_ROIS2) [9] <- "PE-feedback"
names(betas_ROIS2) [10] <- "feed_tau2"

clusters <- kmeans(betas_contrasts_ROIs[,1:3], 5)
clusters <- kmeans(betas_contrasts_ROIS2[,1:3], 5)
betas_contrasts_ROIs$cluster <- as.factor(clusters$cluster)
betas_contrasts_ROIS2$cluster <- as.factor(clusters$cluster)

names(betas_contrasts_ROIs) [7] <- "SITK"
names(betas_contrasts_ROIs) [8] <- "SITR"
names(betas_contrasts_ROIs) [9] <- "SKTR"
names(betas_contrasts_ROIs) [10] <- "SKTK"
names(betas_contrasts_ROIs) [11] <- "IRvsIK"

names(betas_contrasts_ROIS2) [7] <- "SITK"
names(betas_contrasts_ROIS2) [8] <- "SITK_tau2"
names(betas_contrasts_ROIS2) [9] <- "SITR"
names(betas_contrasts_ROIS2) [10] <- "SITR_tau2"
names(betas_contrasts_ROIS2) [11] <- "SKTR"
names(betas_contrasts_ROIS2) [12] <- "SKTR_tau2"
names(betas_contrasts_ROIS2) [13] <- "SKTK"
names(betas_contrasts_ROIS2) [14] <- "SKTK_tau2"
names(betas_contrasts_ROIS2) [15] <- "IRvsIK"
names(betas_contrasts_ROIS2) [16] <- "IRvsIK_tau2"


ggplot(betas_ROIS[,4:6], aes(x = betas_ROIS$X4, y = betas_ROIS$X5, color = betas_ROIS$cluster)) + geom_point() #axial view: #1: left insula/front operculum; #2: right striatum; #3: midbrain/thalamus; #4: left striatum; #5: Left fusiform, lingual, inferior occipital gyrus (V3, V4)
ggplot(betas_ROIS2[,4:6], aes(x = betas_ROIS2$X4, y = betas_ROIS2$X5, color = betas_ROIS2$cluster)) + geom_point() 
ggplot(betas_contrasts_ROIs[,4:6], aes(x = betas_contrasts_ROIs$X4, y = betas_contrasts_ROIs$X5, color = betas_contrasts_ROIs$cluster)) + geom_point() #axial view: #1: Left fusiform, lingual, inferior occipital gyrus (V3, V4); #2: right striatum; #3: left insula/front operculum; #4: midbrain/thalamus; #5: left striatum
ggplot(betas_contrasts_ROIS2[,4:6], aes(x = betas_contrasts_ROIS2$X4, y = betas_contrasts_ROIS2$X5, color = betas_contrasts_ROIS2$cluster)) + geom_point() #axial view: #1: Left fusiform, lingual, inferior occipital gyrus (V3, V4); #2: right striatum; #3: left insula/front operculum; #4: midbrain/thalamus; #5: left striatum


betas_ROIS$clust_label[betas_ROIS$cluster == 1] <- "Linsula_FRNToperc"
betas_ROIS$clust_label[betas_ROIS$cluster == 2] <- "Rstriatum"
betas_ROIS$clust_label[betas_ROIS$cluster == 3] <- "midbrain_thalamus"
betas_ROIS$clust_label[betas_ROIS$cluster == 4] <- "Lstriatum"
betas_ROIS$clust_label[betas_ROIS$cluster == 5] <- "Lfusiform_OCC"
names(betas_ROIS) [7] <- "BETA_PE_DEC"
names(betas_ROIS) [8] <- "BETA_PE_FEED"

betas_ROIS2$clust_label[betas_ROIS2$cluster == 4] <- "Linsula_FRNToperc"
betas_ROIS2$clust_label[betas_ROIS2$cluster == 3] <- "Rstriatum"
betas_ROIS2$clust_label[betas_ROIS2$cluster == 1] <- "midbrain_thalamus"
betas_ROIS2$clust_label[betas_ROIS2$cluster == 2] <- "Lstriatum"
betas_ROIS2$clust_label[betas_ROIS2$cluster == 5] <- "Lfusiform_OCC"

betas_contrasts_ROIs$clust_label[betas_contrasts_ROIs$cluster == 1] <- "Lfusiform_OCC"
betas_contrasts_ROIs$clust_label[betas_contrasts_ROIs$cluster == 2] <- "Rstriatum"
betas_contrasts_ROIs$clust_label[betas_contrasts_ROIs$cluster == 3] <- "Linsula_FRNToperc"
betas_contrasts_ROIs$clust_label[betas_contrasts_ROIs$cluster == 4] <- "midbrain_thalamus"
betas_contrasts_ROIs$clust_label[betas_contrasts_ROIs$cluster == 5] <- "Lstriatum"

betas_contrasts_ROIS2$clust_label[betas_contrasts_ROIS2$cluster == 2] <- "Lfusiform_OCC"
betas_contrasts_ROIS2$clust_label[betas_contrasts_ROIS2$cluster == 4] <- "Rstriatum"
betas_contrasts_ROIS2$clust_label[betas_contrasts_ROIS2$cluster == 3] <- "Linsula_FRNToperc"
betas_contrasts_ROIS2$clust_label[betas_contrasts_ROIS2$cluster == 5] <- "midbrain_thalamus"
betas_contrasts_ROIS2$clust_label[betas_contrasts_ROIS2$cluster == 1] <- "Lstriatum"

ggplot(betas_contrasts_ROIs[,4:6], aes(x = betas_contrasts_ROIs$X4, y = betas_contrasts_ROIs$X5, color = betas_contrasts_ROIs$clust_label)) + geom_point() #axial view: #1: Left fusiform, lingual, inferior occipital gyrus (V3, V4); #2: right striatum; #3: left insula/front operculum; #4: midbrain/thalamus; #5: left striatum

ggplot(betas_contrasts_ROIS[,4:6], aes(x = betas_contrasts_ROIs$X4, y = betas_contrasts_ROIS$X5, color = betas_contrasts_ROIs$clust_label)) + geom_point() #axial view: #1: Left fusiform, lingual, inferior occipital gyrus (V3, V4); #2: right striatum; #3: left insula/front operculum; #4: midbrain/thalamus; #5: left striatum


b <- betas_ROIS
ggplot(b,aes(x = clust_label, y = BETA_PE_DEC)) + geom_boxplot()
ggplot(b,aes(x = clust_label, y = BETA_PE_FEED)) + geom_boxplot()

c <-betas_contrasts_ROIs
ggplot(c,aes(x = clust_label, y = SITK)) + geom_boxplot()
ggplot(c,aes(x = clust_label, y = SITR)) + geom_boxplot()

decPEs <- betas_ROIS[c(1:7,9,10)]
feedPEs <- betas_ROIS[c(1:6,8:10)]
decPEs$align <- "dec"
decPEs$voxels <- 1:1111
names(decPEs) [7] <- "betas"
feedPEs$align <- "feed"
feedPEs$voxels <- 1:1111
names(feedPEs) [7] <- "betas"
PEs <- bind_rows(decPEs, feedPEs, .id = )
names(PEs) [9] <- "ROI"
attach(PEs)
PEs$ROI_rename[PEs$ROI == "Lstriatum"] <- ""
PEs$ROI_rename[PEs$ROI == "Lstriatum"] <- "Left Striatum"
PEs$ROI_rename[PEs$ROI == "Rstriatum"] <- "Right Striatum"
PEs$align_rename <- NA
PEs$align_rename[PEs$align == "dec"] <- "PE-decision"
PEs$align_rename[PEs$align == "feed"] <- "PE-feedback"

ggplot(PEs[PEs$ROI == "Lstriatum" | PEs$ROI == "Rstriatum",],aes(x = align_rename, y = betas)) + geom_jitter() + facet_wrap(~ROI_rename)+
  scale_color_manual(values = c("black", "gray60"))+
  theme_bw()+
  xlab("Alignment of PEs") + ylab("Estimated Mean of Beta Coefficients")
ggsave("pe_jitter_plots.pdf", width = 4.5, height = 3)


m1 <- lme4::lmer(data = PEs[PEs$ROI == "Lstriatum" | PEs$ROI == "Rstriatum",], betas ~ align+ROI+(1|voxels), na.action = na.omit)
#Plots of betas + clusters
lsm <- emmeans::lsmeans(m1,~align|ROI)
em <- cld(lsm, alpha = 0.05, Letters = letters, adjust = "tukey", interaction = "TRUE")
em$ROI = revalue(em$ROI,c("Lstriatum" = "Left Striatum", "Rstriatum" = "Right Striatum"))
em$align = revalue(em$align,c("dec" = "PE-Decision", "feed" = "PE-Feedback"))
pd = position_dodge(0.5)
ggplot(data = em, aes(x=align, y=lsmean, color = ROI)) +
  geom_errorbar(aes(ymin = lsmean-SE, ymax = lsmean+SE), position = pd, width = .1) +
  geom_point(position=pd, size = 2)+
  scale_color_manual(values = c("black", "gray60"))+
  theme_bw()+
  xlab("Alignment of PEs") + ylab("Estimated Mean of Beta Coefficients")
ggsave("pes_at_alignment.pdf", width = 4.5, height = 3)

b <- betas_ROIS2
ggplot(b,aes(x = clust_label, y = b$`PE-decision`)) + geom_boxplot()
ggplot(b,aes(x = clust_label, y = b$`PE-feedback`)) + geom_boxplot()

M_PEdec_LS <- mean(b$`PE-decision`[b$clust_label=="Lstriatum"])
T_PEdec_LS <- mean(sqrt(b$dec_tau2[b$clust_label=="Lstriatum"]))
M_PEdec_RS <- mean(b$`PE-decision`[b$clust_label=="Rstriatum"])
T_PEdec_RS <- mean(sqrt(b$dec_tau2[b$clust_label=="Rstriatum"]))
M_PEfeed_LS <- mean(b$`PE-feedback`[b$clust_label=="Lstriatum"])
T_PEfeed_LS <- mean(sqrt(b$feed_tau2[b$clust_label=="Lstriatum"]))
M_PEfeed_RS <- mean(b$`PE-feedback`[b$clust_label=="Rstriatum"])
T_PEfeed_RS <- mean(sqrt(b$feed_tau2[b$clust_label=="Rstriatum"]))
SD_PE_LS <- mean(cbind(sqrt(b$dec_tau2[b$clust_label=="Lstriatum"]), sqrt(b$feed_tau2[b$clust_label=="Lstriatum"])))
SD_PE_RS <- mean(cbind(sqrt(b$dec_tau2[b$clust_label=="Rstriatum"]), sqrt(b$feed_tau2[b$clust_label=="Rstriatum"])))

left <- z.test(b$`PE-decision`[b$clust_label=="Lstriatum"], b$`PE-feedback`[b$clust_label=="Lstriatum"], mu = 0, sigma.x = T_PEdec_LS, sigma.y = T_PEfeed_LS, conf.level = 0.95)
right <- z.test(b$`PE-decision`[b$clust_label=="Rstriatum"], b$`PE-feedback`[b$clust_label=="Rstriatum"], mu = 0, sigma.x = T_PEdec_RS, sigma.y = T_PEfeed_RS, conf.level = 0.95)

PEdecfeed.data <- data.frame(
  align = c ("PE-Decision", "PE-Feedback", "PE-Decision", "PE-Feedback"), 
  ROI = c("Left Striatum","Left Striatum","Right Striatum","Right Striatum"),
  mean = c(M_PEdec_LS,M_PEfeed_LS,M_PEdec_RS,M_PEfeed_RS), 
  tau = c(T_PEdec_LS, T_PEfeed_LS, T_PEdec_RS, T_PEfeed_RS),
  SEtau = c(T_PEdec_LS/sqrt(22), T_PEfeed_LS/sqrt(22), T_PEdec_RS/sqrt(22), T_PEfeed_RS/sqrt(22))
)

ggplot(data = PEdecfeed.data, aes(x=align, y=mean)) +
  geom_errorbar(aes(ymin = mean-SEtau, ymax = mean+SEtau), position = pd, width = .1) +
  geom_point(position=pd, size = 2)+
  scale_color_manual(values = c("black", "gray60"))+
  theme_bw()+
  facet_wrap(~ROI)+
  xlab("Alignment of PEs") + ylab("Mean of Beta Coefficients")
ggsave("pes_at_alignment2.pdf", width = 4.5, height = 4)


c <- betas_contrasts_ROIS2
M_SITK_LS <- mean(c$SITK[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"])
T_SITK_LS <- mean(sqrt(c$SITK_tau2[c$clust_label=="Lstriatum"]))
M_SITR_LS <- mean(c$SITR[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"])
T_SITR_LS <- mean(sqrt(c$SITR_tau2[c$clust_label=="Lstriatum"]))
M_SKTR_LS <- mean(c$SKTR[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"])
T_SKTR_LS <- mean(sqrt(c$SKTR_tau2[c$clust_label=="Lstriatum"]))
M_SKTK_LS <- mean(c$SKTK[c$clust_label=="Lstriatum"])
T_SKTK_LS <- mean(sqrt(c$SKTR_tau2[c$clust_label=="Lstriatum"]))

M_SITK_RS <- mean(c$SITK[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"])
T_SITK_RS <- mean(sqrt(c$SITK_tau2[c$clust_label=="Rstriatum"]))
M_SITR_RS <- mean(c$SITR[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"])
T_SITR_RS <- mean(sqrt(c$SITR_tau2[c$clust_label=="Rstriatum"]))
M_SKTR_RS <- mean(c$SKTR[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"])
T_SKTR_RS <- mean(sqrt(c$SKTR_tau2[c$clust_label=="Rstriatum"]))
M_SKTK_RS <- mean(c$SKTK[c$clust_label=="Rstriatum"])
T_SKTK_RS <- mean(sqrt(c$SKTR_tau2[c$clust_label=="Rstriatum"]))

OutcomeFeedback.data <- data.frame(
  outcometype = c ("S-Invested\nT-Kept", "S-Invested\nT-Returned", "S-Kept\nT-Kept", "S-Kept\nT-Returned",
                   "S-Invested\nT-Kept", "S-Invested\nT-Returned", "S-Kept\nT-Kept", "S-Kept\nT-Returned"), 
  ROI = c("Left Striatum","Left Striatum","Left Striatum","Left Striatum","Right Striatum","Right Striatum","Right Striatum","Right Striatum"),
  mean = c(M_SITK_LS,M_SITR_LS,M_SKTK_LS,M_SKTR_LS,M_SITK_RS,M_SITR_RS,M_SKTK_RS,M_SKTR_RS), 
  tau = c(T_SITK_LS,T_SITR_LS,T_SKTK_LS,T_SKTR_LS,T_SITK_RS,T_SITR_RS,T_SKTK_RS,T_SKTR_RS),
  SEtau = c(T_SITK_LS/sqrt(22),T_SITR_LS/sqrt(22),T_SKTK_LS/sqrt(22),T_SKTR_LS/sqrt(22),T_SITK_RS/sqrt(22),T_SITR_RS/sqrt(22),T_SKTK_RS/sqrt(22),T_SKTR_RS/sqrt(22))
)

ggplot(data = OutcomeFeedback.data, aes(x=outcometype, y=mean)) +
  geom_errorbar(aes(ymin = mean-SEtau, ymax = mean+SEtau), position = pd, width = .1) +
  geom_point(position=pd, size = 2)+
  scale_color_manual(values = c("black", "gray60"))+
  theme_bw()+
  facet_wrap(~ROI)+
  xlab("Outcome Type") + ylab("Mean of Beta Coefficients")
ggsave("outcometype2.pdf", width = 6.5, height = 4)

SITKvsSKTK_LS <- z.test(c$SITK[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], c$SKTK[c$clust_label=="Lstriatum"], mu = 0, sigma.x = T_SITK_LS, sigma.y = T_SKTK_LS, conf.level = 0.95)
SITKvsSITR_LS <- z.test(c$SITK[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], c$SITR[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], mu = 0, sigma.x = T_SITK_LS, sigma.y = T_SITR_LS, conf.level = 0.95)
SITKvsSKTR_LS <- z.test(c$SITK[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], c$SKTR[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], mu = 0, sigma.x = T_SITK_LS, sigma.y = T_SKTR_LS, conf.level = 0.95)
SITRvsSKTK_LS <- z.test(c$SITR[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], c$SKTK[c$clust_label=="Lstriatum"], mu = 0, sigma.x = T_SITR_LS, sigma.y = T_SKTK_LS, conf.level = 0.95)
SKTRvsSKTK_LS <- z.test(c$SKTR[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], c$SKTK[c$clust_label=="Lstriatum"], mu = 0, sigma.x = T_SKTR_LS, sigma.y = T_SKTK_LS, conf.level = 0.95)
SITRvsSKTR_LS <- z.test(c$SITR[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], c$SKTR[c$clust_label=="Lstriatum"]+c$SKTK[c$clust_label=="Lstriatum"], mu = 0, sigma.x = T_SITR_LS, sigma.y = T_SKTR_LS, conf.level = 0.95)


SITKvsSKTK_RS <- z.test(c$SITK[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], c$SKTK[c$clust_label=="Rstriatum"], mu = 0, sigma.x = T_SITK_RS, sigma.y = T_SKTK_RS, conf.level = 0.95)
SITKvsSITR_RS <- z.test(c$SITK[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], c$SITR[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], mu = 0, sigma.x = T_SITK_RS, sigma.y = T_SITR_RS, conf.level = 0.95)
SITKvsSKTR_RS <- z.test(c$SITK[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], c$SKTR[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], mu = 0, sigma.x = T_SITK_RS, sigma.y = T_SKTR_RS, conf.level = 0.95)
SITRvsSKTK_RS <- z.test(c$SITR[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], c$SKTK[c$clust_label=="Rstriatum"], mu = 0, sigma.x = T_SITR_RS, sigma.y = T_SKTK_RS, conf.level = 0.95)
SKTRvsSKTK_RS <- z.test(c$SKTR[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], c$SKTK[c$clust_label=="Rstriatum"], mu = 0, sigma.x = T_SKTR_RS, sigma.y = T_SKTK_RS, conf.level = 0.95)
SITRvsSKTR_RS <- z.test(c$SITR[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], c$SKTR[c$clust_label=="Rstriatum"]+c$SKTK[c$clust_label=="Rstriatum"], mu = 0, sigma.x = T_SITR_RS, sigma.y = T_SKTR_RS, conf.level = 0.95)
