rm(list = ls())

############Arguments
#Model:
Model_PEValue<-T
Model_PEValue_CusCmat<-F
Model_PEonly<-F
###############
#Sub-set
run_on_noComputer<-F
run_on_onlyComputer<-F
run_on_HConly<-F
###############
run_step=3:4
lvl3_covariate_names = c("Intercept","Group_LETH")
bad_ids = c(222027, 440006, 221846, 221945)

require("devtools")
devtools::install_github("DecisionNeurosciencePsychopathology/fMRI_R",force = F)
devtools::install_github("PennStateDEPENdLab/dependlab",force=F)
library(fslpipe)
######Functions: #####
require("dplyr")
prep_trust_fsl<-function(dfx=NULL){
  if(is.null(dfx)){stop("NO DATA")}
  dfx<-dfx[order(as.numeric(dfx$TrialNumber)),]
  dfx$dec_bin<-ifelse(dfx$PartDecides=="share",1,-1)
  dfx$trustdec_bin<-ifelse(dfx$TrusteeDecides=="share",1,-1)
  dfx$QC<-dfx$dec_bin
  dfx$congruent<-ifelse(dfx$PartDecides == dfx$TrusteeDecides,1,-1)
  dfx$bothShare<-ifelse(dfx$PartDecides=="share" & dfx$TrusteeDecides=="share",1,-1)

  dfx$trustee_bin<-NA
  dfx$trustee_bin[which(dfx$identity=="neutral")]<- 0
  dfx$trustee_bin[which(dfx$identity=="bad")]<- -1
  dfx$trustee_bin[which(dfx$identity=="good")]<- 1
  #Alex says temporarily we set computer to be bad;
  dfx$trustee_bin[which(dfx$identity=="computer")]<- -1

  #Trustee
  dfx$trustee_Good<-0;dfx$trustee_Bad<-0;dfx$trustee_Neutral<-0;dfx$trustee_PC<-0;
  dfx$trustee_Good[which(dfx$identity=="good")]<- 1
  dfx$trustee_Bad[which(dfx$identity=="bad")]<- 1
  dfx$trustee_Neutral[which(dfx$identity=="neutral")]<- 1
  dfx$trustee_PC[which(dfx$identity=="computer")]<- 1

  dfx$PE<-dplyr::lead(dfx$PE_1)
  dfx$V_post<-dplyr::lead(dfx$V_1)
  dfx$V_pre <- dfx$V_1

  ##Policy Model:
  dfx$PE_policy<-dfx$PE
  dfx$PE_policy[which(dfx$PartDecides=="keep")]<- (-dfx$PE[which(dfx$PartDecides=="keep")])

  dfx$V_policy<-dfx$V_post
  dfx$V_policy[which(dfx$PartDecides=="keep")]<- (-dfx$V_post[which(dfx$PartDecides=="keep")])


  # dfx$PChosen_post<-dfx$V_postupdate
  # dfx$PChosen_post[which(dfx$PartDecides=="keep")]<-1-(dfx$V_postupdate[which(dfx$PartDecides=="keep")])
  #
  # dfx$PChosen_pre<-dfx$V_1
  # dfx$PChosen_pre[which(dfx$PartDecides=="keep")]<-1-(dfx$PChosen_pre[which(dfx$PartDecides=="keep")])
  #
  # dfx$VShare_utility_post<-(dfx$V_postupdate * 1.5) -1
  # dfx$VShare_utility_pre<-(dfx$V_1 * 1.5) -1
  #
  # dfx$PE_chosen<-dfx$PE_1
  # dfx$PE_chosen[which(dfx$PartDecides=="keep")]<- (-dfx$PE_1[which(dfx$PartDecides=="keep")])
  #
  # dfx$V_Chosen_pre<-dfx$V_1
  # dfx$V_Chosen_pre[which(dfx$PartDecides=="keep")]<- -(dfx$V_1[which(dfx$PartDecides=="keep")])

  dfx$computer<-0
  dfx$computer[which(dfx$identity=="computer")]<-1

  dfx$Run<-1
  dfx$Trial<-as.numeric(dfx$TrialNumber)

  dfx_sp<-split(dfx,dfx$identity)

  trust_df<-do.call(rbind,lapply(dfx_sp,function(rjx){
    rjxa<-data.frame(onset=(as.numeric(rjx$partnerchoice.OnsetTime)[1] - as.numeric(rjx$firstFixation.OnsetTime[1]) )/1000,
                     duration=(as.numeric(rjx$ITIfixation.OffsetTime)[nrow(rjx)] - as.numeric(rjx$partnerchoice.OnsetTime)[1])/1000,
                     run=1)
    rjxa$computer<-unique(rjx$computer)
    rjxa$trustee_bin<-unique(rjx$trustee_bin)
    return(rjxa)
  }))
  trust_df$trial<-1:nrow(trust_df)

  #Set up events
  evt_list<-list(Decision=data.frame(event="Decision",
                                     onset= ( as.numeric(dfx$partnerchoice.OnsetTime) - as.numeric(dfx$firstFixation.OnsetTime[1]) )/1000  ,
                                     duration=as.numeric(dfx$partnerchoice.RT)/1000,
                                     run=dfx$Run,
                                     trial=dfx$Trial),
                 Outcome=data.frame(event="Outcome",
                                    onset=(as.numeric(dfx$outcome.OnsetTime) - as.numeric(dfx$firstFixation.OnsetTime[1]) )/1000,
                                    duration=as.numeric(dfx$outcome.OffsetTime)/1000 - as.numeric(dfx$outcome.OnsetTime)/1000,
                                    run=dfx$Run,
                                    trial=dfx$Trial),
                 # TrustBlock=data.frame(event="TrustBlock",
                 #                       onset=trust_df$onset,
                 #                       duration=trust_df$duration,
                 #                       run=trust_df$run,
                 #                       trial=trust_df$trial),
                 QC=data.frame(event="QC",
                               onset=(as.numeric(dfx$partnerchoice.OnsetTime) - as.numeric(dfx$firstFixation.OnsetTime[1]) )/1000,
                               duration=as.numeric(dfx$partnerchoice.RT)/1000,
                               run=dfx$Run,
                               trial=dfx$Trial)
  )
  for (i in 1:length(evt_list)) {
    if (i==1) {ktz<-evt_list[[i]]} else {
      ktz<-rbind(ktz,evt_list[[i]])}
  }
  evt_list[["allconcat"]]<-ktz

  value<-as.list(dfx)
  #value$Trustee_block<-trust_df$trustee_bin
  #value$Computer_block<-trust_df$computer
  output<-list(event.list=evt_list,output.df=dfx,value=value)
  return(output)

}

#Setting up FSL global enviroment variables in case we are using RStudio
fsl_2_sys_env()
load("trust_df.rdata")

######
#Actual arguments for each model. Should follow template: github.com/DecisionNeurosciencePsychopathology/fMRI_R
####BE AWARE!
argu<-as.environment(list(nprocess=12,run_steps=run_step,forcereg=F,
                          cfgpath="/Volumes/bek/autopreprocessing_pipeline/Learn/trust.cfg",
                         func.nii.name="nfswudktm*[0-9]_[0-9].nii.gz",
                         regtype=".1D", adaptive_gfeat=TRUE,adaptive_ssfeat=TRUE,centerscaleall=FALSE,
                         templatedir="/usr/local/ni_tools/standard_templates/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_brain_2.3mm.nii",
                         ####Output Paths wherever you want them to be;
                         regpath="/gpfs/group/LiberalArts/default/mnh5174_collab/bsocial/trust/reg",
                         subj_outputroot="/gpfs/group/LiberalArts/default/mnh5174_collab/bsocial/trust/ssub_analysis",
                         glvl_output="/gpfs/group/LiberalArts/default/mnh5174_collab/bsocial/trust/grp_analysis"
))
#DO NOT PROVIDE THOSE TWO AND IT WILL BE FINE;
#We are refiting lvl2
argu$lvl2_overwrite<-F
argu$lvl3_overwrite<-F
argu$lvl1_retry<-F
#argu$lvl3_type<-"flame"
#argu$thresh_cluster_extent<-3.1
#argu$thresh_cluster_mass<-3.1
argu$cfg<-cfg_info(cfgpath = argu$cfgpath)
argu$ss_zthreshold<-1.96  #This controls the single subject z threshold (if enabled in template)
argu$ss_pthreshold<-0.05 #This controls the single subject p threshold (if enabled in template)

#PBS argu
argu$run_on_pbs<-FALSE

###############
###Data processing;
trust_sp<-split(trust_merged,trust_merged$ID)
trust_sp<-trust_sp[!names(trust_sp) %in% bad_ids]

if(run_on_noComputer){
  hasComputer_subj<-sapply(trust_sp,function(x){unique(x$hasComputer)})
  trust_sp<-trust_sp[!hasComputer_subj]
  message("Run on subjects who don't have computer runs...The N is ",length(trust_sp))
}
if(run_on_onlyComputer){
  hasComputer_subj<-sapply(trust_sp,function(x){unique(x$hasComputer)})
  trust_sp<-trust_sp[hasComputer_subj]
  message("Run on subjects who  have computer runs...The N is ",length(trust_sp))
}
if(run_on_HConly){
  HCsubj<-sapply(trust_sp,function(x){"HC" %in% x$Group})
  trust_sp<-trust_sp[HCsubj]
  argu$run_these_ID <- names(trust_sp)
  message("Run on subjects who are healthy controls...The N is ",length(trust_sp))
}


trust_group_df<-do.call(rbind,lapply(trust_sp, function(x){unique(x[c("ID","Group","Lethality")])}))
trust_group_df$Group_LETH<-ifelse(trust_group_df$Group=="ATT",toupper(trust_group_df$Lethality),trust_group_df$Group)
trust_group_df$Group_LETH<-factor(trust_group_df$Group_LETH,levels = c("HL","LL","HC","NON"))
trust_group_df$Group<-factor(trust_group_df$Group,levels = c("ATT","HC","NON"))
trust_sp_f<-lapply(trust_sp, function(x){list(dfx=x)})

if(Model_PEValue) {
  argu$model_name="PEValue"
  argu$gridpath="./grid/grid_PEValue.csv"
  argu$centerscaleall=TRUE
  argu$lvl3_ref_df<-trust_group_df
  argu$lvl3_covariate_names<-lvl3_covariate_names
}

if(Model_PEonly) {
  argu$model_name="PE_only"
  argu$gridpath="./grid/grid_PE.csv"
  argu$centerscaleall=TRUE
  argu$lvl3_ref_df<-trust_group_df
  argu$lvl3_covariate_names<-lvl3_covariate_names
}

if(Model_PEValue_CusCmat) {
  argu$model_name="PEValue_Int"
  if(run_on_noComputer) {
    argu$gridpath="./grid/grid_PEValue_Interaction_nopc.csv"
  }else if (run_on_onlyComputer) {
    argu$gridpath="./grid/grid_PEValue_Interaction_pconly.csv"
    }

  argu$centerscaleall=TRUE
  argu$lvl3_ref_df<-trust_group_df
  argu$lvl3_covariate_names<-lvl3_covariate_names


  argu$dsgrid<-read.csv(argu$gridpath,stringsAsFactors = F)
  IndexDF<-data.frame(evanme=argu$dsgrid$name,MainContrast=sapply(strsplit(argu$dsgrid$name,"_"),`[[`,2),Term=sapply(strsplit(argu$dsgrid$name,"_"),`[[`,1),stringsAsFactors = F)
  TrusteeNum<-3
  #First time using the customized contrast
  argu$lvl1_cmat<-do.call(rbind,lapply(unique(IndexDF$MainContrast),function(mainx){
    mxa<-matrix(data = 0,nrow = (nrow(IndexDF[IndexDF$MainContrast==mainx,])+1),ncol = nrow(IndexDF))
    colnames(mxa)<-IndexDF$evanme
    rownumx<-which(IndexDF$MainContrast==mainx)
    mxa[1:length(rownumx),rownumx]<-diag(x = 1,nrow = length(rownumx),ncol = length(rownumx))
    mxa[length(rownumx)+1,rownumx] <- 1 / length(rownumx)
    rownames(mxa)<-c(IndexDF$evanme[rownumx],paste("Overall",mainx,sep = "_"))

    exgr<-t(combn(x = IndexDF[rownumx,"Term"],m = 2))
    mxb<-matrix(data = 0,nrow = nrow(exgr),ncol = nrow(IndexDF))
    colnames(mxb)<-IndexDF$evanme
    rownames(mxb)<-paste(apply(exgr,1,paste,collapse="_grt_"),mainx,sep = "_")

    for (numx in 1:nrow(exgr)) {
      x<-exgr[numx,]
      mxb[numx,which(IndexDF$MainContrast==mainx & IndexDF$Term==x[1])]<-1
      mxb[numx,which(IndexDF$MainContrast==mainx & IndexDF$Term==x[2])]<-(-1)
    }
    mxc<-rbind(mxa,mxb)
    return(mxc)
  }))
  argu$lvl1_cmat<-argu$lvl1_cmat[which(!grepl("TPC",rownames(argu$lvl1_cmat))),]
}


argu$exclude_these_ID<-bad_ids

#trust_sp_f<-trust_sp_f[1]
#trust_sp<-trust_sp[c(2,3)]
fsl_pipe(
  argu=argu, #This is the arguments environment, each model should have a different one;
  prep.call.func="prep_trust_fsl", #This should be a character string that's the name of the prep proc function
  prep.call.allsub=trust_sp_f #List of ID list of arguments for prep.call.
)

