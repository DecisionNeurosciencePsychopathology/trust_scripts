trust_vba_df<-read.csv("~/Box/skinner/projects_analyses/Trust/mfx_analyses/learn_policy_FixedPara_no_multisession_FourTrusteekappaST_mfx_trial_outputs_by_timestep.csv",stringsAsFactors = F)

library(rprime)
rootdir<-"~/Box/skinner/data/eprime/scan_trust"
#learn<-new.env()
#load("~/Box/skinner/data/RedCap Data/BSocial/Redcap.bsocial.rdata",envir = bsocial)
# loaded manually
availableIDs<-list.dirs(rootdir,recursive = F,full.names = F)

trust_vba_sp <- split(trust_vba_df,trust_vba_df$id)

trust_behav_ls<-lapply(trust_vba_sp, function(gx) {
  IDx<-unique(gx$id)
  print(IDx)
  if(!IDx %in% availableIDs) {message(IDx," has no behavioral data.");return(NULL)}

  rx<-list.files(file.path(rootdir,IDx),pattern = "trust.*.mat",full.names = T)
  if(length(rx)!=1) {message(IDx," has no behavioral data. No .mat file");return(NULL)}

  ga<-R.matlab::readMat(rx)
  if(is.null(ga$b)){message(IDx," has no behavioral data. Mis-config mat file");return(NULL)}
  ga <- ga$b
  ga<-lapply(drop(ga),function(x){
    x <- lapply(drop(x),drop)
    if( length(x) == length(unlist(x)) ) {
      x <- unlist(x,recursive = T,use.names = F)
    }
    return(x)
  })



  ga[paste("trustee",1:4,sep = "")]<-NULL

  ga<-lapply(ga,function(xa){
    if(!any(!sapply(xa,is.null))){return(NA)}
    if(length(xa)!=192 && length(xa)!=1) {xa[(length(xa)+1):192]<-NA}
    return(xa)
  })

  ga_df<-as.data.frame(ga,stringsAsFactors = F)
  if(any(is.na(ga_df$firstFixation.OnsetTime))){ga_df$firstFixation.OnsetTime <- (ga_df$partnerchoice.OnsetTime[1]-250)}
  ga_df$Run <-match(ga_df$identity,unique(ga_df$identity))
  ga_df$Run[is.na(ga_df$identity)]<-NA
  if(nrow(ga_df)!=192){message(IDx," has behavioral data that has ", nrow(ga_df)," rows")}
  if(nrow(gx)!=192){message(IDx," has vba data that has ", nrow(gx)," rows")}
  #Goal here is to get the correct onset time for block 2,3,4
  ga_index<-which(!duplicated(ga_df$identity))
  ga_index<- ga_index[ga_index!=1]
  ga_df$firstFixation.OnsetTime [ga_df$Run!=1 & !is.na(ga_df$Run)] <- ga_df$ITIfixation.OnsetTime[na.omit((ga_index-1)[ga_df$Run[ga_df$Run!=1]-1  ])]


  ga_mg<-merge(ga_df,gx,by.x = "TrialNumber",by.y = "asc_trial",all = T)
  ga_mg$ID<-IDx
  ga_mg[ga_mg==-999] <- NA



  if(any(na.omit(as.numeric(ga_mg$PartDecides=="share") != ga_mg$y_1))){message("This person's data does not match: ",unique(ga_mg$ID));return(NULL)}

  return(ga_mg)
})

allnames<-unique(unlist(unique(lapply(trust_behav_ls,names))))
trust_behav_ls<-lapply(trust_behav_ls,function(xa){xa[allnames[!allnames %in% names(xa)]]<-NA;return(xa)})

trust_merged_a<-do.call(rbind,trust_behav_ls)

masterdemo <- bsrc::bsrc.checkdatabase2(protocol = ptcs$masterdemo)
idmap<-masterdemo$data[c("registration_redcapid","registration_wpicid","registration_edu","registration_gender","registration_dob","registration_group","registration_lethality")]
names(idmap) <- c("masterdemo_id","wpic_id","edu","gender","dob","Group","Lethality")
trust_merged<-bsrc.findid(trust_merged_a,idmap = idmap,id.var = "ID")
save(trust_merged,file = "./trust_cluster/trust_learn_df.rdata")
paste0("rsync -a -P ",getwd(),"/trust_cluster/. ayd5415@aci-b.aci.ics.psu.edu:~/trust_cluster/.")


allinfo<-do.call(rbind,lapply(list.dirs("~/learn_proc/MR_Proc/",recursive = F),function(xa){
  daz<-do.call(cbind,lapply(1:4,function(n){
    da<-data.frame()
    da[1,paste0("trust",n)]<-file.exists(file.path(xa,"mni_7mm_aroma",paste0("trust",n),paste0("nfaswuktm_trust",n,"_7.nii.gz")))
    return(da)
  }))
  daz$ID <- basename(xa)
  return(daz)
}))



firstlvl<-c(202200,       210701,       214082,       218219,       219886,       220186,       220597,       220927,       221181,       35780,
  202278,       211038,       215088,       218572,       219944,       220244,       220678,       220947,       221183,       46069,
  203264,       211101,       215537,       218714,       219949,       220299,       220684,       220963,       221206,       881075,
  203803,       211253,       215622,       219044,       219954,       220399,       220740,       220976,       221212,       881105,
  204015,       212019,       216011,       219392,       220043,       220477,       220758,       221018,       221226,       881132,
  206270,       212020,       217008,       219619,       220051,       220513,       220789,       221036,       221261,       881224,
  207224,       212794,       217048,       219658,       220064,       220523,       220865,       221050,       221273,       881230,
  208510,       213163,       217173,       219676,       220104,       220531,       220889,       221099,       221292,
210290,       213704,       217293,       219772,       220152,       220553,       220913,       221140,       221295,
  210374,       213868,       217988,       219809,       220184,       220566,       220926,       221164,       221298,
220989,221085
)
unique(trust_merged$ID[!trust_merged$ID %in% firstlvl])


trust_sp<-split(trust_merged,trust_merged$ID)

any(sapply(trust_sp,function(b){
  out<-prep_trust_fsl(b)
  any(na.omit(out$event.list$Outcome$duration)<0)
}))



p_t<-c(1.802,3.72,4.187,3.085)*1000



for ( block in 1:4){
  message("####",block)
  start_time <- b_sp[[block]]$outcome.OnsetTime[1] - p_t[block]
  message("Starting at:", start_time)
  print(b_sp[[block]] [1,grepl("onset",tolower(names(b_sp[[block]])))] - start_time)
  message("cal_matlab:" )
}



emptyGX<-is.na(b$ITIfixation.OnsetTime) | b$ITIfixation.OnsetTime==(-999)

if (any(emptyGX)) {
  message("A")
  b$ITIfixation_OnsetTime[emptyGX]=-999;
  b$ITIfixation_OffsetTime[145:192]= -999;
}

if (b$ITIfixation_OnsetTime[192] == -999 ) {
  message("B if")
  firstfix_Onset = b$firstFixation_OnsetTime[1];
  trial2_ITI = 1;
  trial48_ITI = 47;
} else {
  message("B else")
  firstfix_Onset = b$ITIfixation_OnsetTime[1];
  trial2_ITI = 2;
  trial48_ITI = 48;
}

if iscell(b.partnerchoice_OffsetTime)
  partnerchoice_OffsetTime=b.displaychoice_OnsetTime-1;
else
  partnerchoice_OffsetTime=b.partnerchoice_OffsetTime;
end

fixations = [];
fixations.event_beg=zeros(48,4);
fixations.event_end=zeros(48,4);
taskness.event_beg =zeros(48,4);
taskness.event_end=zeros(48,4);
decision.event_beg=zeros(48,4);
decision.event_end=zeros(48,4);
feedback.event_beg=zeros(48,4);
feedback.event_end=zeros(48,4);
if iscell(b.partnerchoice_RESP)
b.missed_trials = b.decisions == 0;
else
  b.missed_trials = (b.partnerchoice_RESP==-999);
end
b.notmissed_trials = ~(b.missed_trials);
b.keepNmiss_trials = (b.decisions==-1 | b.missed_trials);
b.shareNmiss_trials = (b.decisions==1 | b.missed_trials);


trial1_index = 1;
trial48_index = 48;



#
# loading_ps<-parallel::makeCluster(parallel::detectCores(),outfile="",type = "FORK")
# drx<-parallel::parSapply(loading_ps,list.dirs(rootdir,recursive = F,full.names = F), function(dirj) {
#   rx<-list.files(file.path(rootdir,dirj),pattern = "^trust.*_scan_.*.txt$",full.names = T)
#   if(length(rx)==1){
#     test1<-rprime::read_eprime(rx)
#     test2<-rprime::extract_chunks(test1) #chuck it
#     tdata<-FrameList(test1) #tdata now is a list
#     #preview_levels(tdata)
#     onlytwo<-keep_levels(tdata,3)
#     testdf<-to_data_frame(onlytwo)
#     if(nrow(testdf)>0){
#       testdf$ID<-dirj
#       #testdf$Group<-bsocial$data$registration_group[match(dirj,bsocial$data$registration_redcapid)]
#       return(testdf)}else{return(NULL)}
#   }else{NULL}
# })
# parallel::stopCluster(loading_ps)
# drx<-drx[!sapply(drx,is.null)]
# save("drx",file = "./trust_cluster/trust_learn_raw.rdata")
# trust_behav_ls<-lapply(drx,function(dra){
#   print(unique(dra$ID))
#   dra$hasComputer<-unique("computer" %in% dra$identity)
#   if(!unique(dra$ID) %in% trust_vba_df$id ){
#     message("This person: ",unique(dra$ID),"does not have a vba data...Skiping...")
#     return(NULL)
#   }
#   dra_v<-trust_vba_df[trust_vba_df$id %in% unique(dra$ID),-1]
#   if(nrow(dra)>192){dra<-dra[1:192,]}
#   if(nrow(dra_v)<192) {
#     message("This person: ",unique(dra$ID),"does not enough trials in their VBA data...adding one to the last row...")
#     dra_v[(nrow(dra_v)+1):192,]<-NA
#   }
#
#   if(nrow(dra_v)>nrow(dra)) {
#     message("This person: ",unique(dra$ID),"does not enough trials in their behavioral dataframe...adding fillers...")
#     dra[(nrow(dra)+1):nrow(dra_v),]<-NA
#   }
#
#   drb<-cbind(dra,dra_v)
#
#   gx<-rep_len(c(1,2,3,4),192)
#   drb$Run<-gx[order(gx)]
#   return(drb)
# })
# allnames<-unique(unlist(unique(lapply(trust_behav_ls,names))))
# trust_behav_ls<-lapply(trust_behav_ls,function(xa){xa[allnames[!allnames %in% names(xa)]]<-NA;return(xa)})
#
# trust_behav_df<-do.call(rbind,trust_behav_ls)

#Getting clinical data:






