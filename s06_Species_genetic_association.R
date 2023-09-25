### Species genetic association
### 2023-09-19
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")

## 1 Preparation
### 1.1 Import
##install.packages("tmvnsim")
source("functions.R")

### 1.2 Inputs
France_info<- read.table("01.cleanData/SV_info/France_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)
USA_UK_info<- read.table("01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

France_basic<-read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")
USA_UK_basic<-read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv")

all_basic<-rbind(France_basic,USA_UK_basic)

PRJEB22863_basic<-subset(all_basic,dataset=="PRJEB22863",drop=T)
PRJEB22863_NSCLC_basic<-subset(PRJEB22863_basic,cancer_type=="NSCLC",drop=T)
PRJEB22863_RCC_basic<-subset(PRJEB22863_basic,cancer_type=="RCC",drop=T)
PRJNA397906_basic<-subset(all_basic,dataset=="PRJNA397906",drop=T)
PRJNA751792_basic<-subset(all_basic,dataset=="PRJNA751792",drop=T)
PRJNA762360_basic<-subset(all_basic,dataset=="PRJNA762360",drop=T)
PRJNA770295_basic<-subset(all_basic,dataset=="PRJNA770295",drop=T)
PRJEB43119_basic<-subset(all_basic,dataset=="PRJEB43119",drop=T)
PRJNA541981_basic<-subset(all_basic,dataset=="PRJNA541981",drop=T)

###  os<-subset(all_basic,is.na(os)==F,drop=T)
###  pfs<-subset(all_basic,is.na(pfs)==F,drop=T)
###  table(os$study)
###  table(pfs$study)
###  response<-subset(all_basic,is.na(response_code)==F,drop=T)
###  table(response$study)

PRJEB22863_abun<-read.table("01.cleanData/mbio_all/PRJEB22863_SV_species_abun.tsv", check.names = F) 
PRJNA397906_abun<-read.table("01.cleanData/mbio_all/PRJNA397906_SV_species_abun.tsv", check.names = F)  
PRJNA751792_abun<-read.table("01.cleanData/mbio_all/PRJNA751792_SV_species_abun.tsv", check.names = F) 
PRJNA762360_abun<-read.table("01.cleanData/mbio_all/PRJNA762360_SV_species_abun.tsv", check.names = F) 
PRJNA770295_abun<-read.table("01.cleanData/mbio_all/PRJNA770295_SV_species_abun.tsv", check.names = F)
PRJEB43119_abun<-read.table("01.cleanData/mbio_all/PRJEB43119_SV_species_abun.tsv", check.names = F)
PRJNA541981_abun<-read.table("01.cleanData/mbio_all/PRJNA541981_SV_species_abun.tsv", check.names = F)

PRJEB22863_NSCLC_abun<-PRJEB22863_abun[PRJEB22863_NSCLC_basic$id,]
PRJEB22863_RCC_abun<-PRJEB22863_abun[PRJEB22863_RCC_basic$id,]

load("01.cleanData/SV_all/distMat/PRJNA397906_msv_dist_std.RData")
PRJNA397906_msv_dist_std<-msv_dist_std
load("01.cleanData/SV_all/distMat/PRJNA762360_msv_dist_std.RData")
PRJNA762360_msv_dist_std<-msv_dist_std
load("01.cleanData/SV_all/distMat/PRJNA770295_msv_dist_std.RData")
PRJNA770295_msv_dist_std<-msv_dist_std
load("01.cleanData/SV_all/distMat/PRJEB43119_msv_dist_std.RData")
PRJEB43119_msv_dist_std<-msv_dist_std
load("01.cleanData/SV_all/distMat/PRJEB22863_msv_dist_std.RData")
PRJEB22863_msv_dist_std<-msv_dist_std
load("01.cleanData/SV_all/distMat/PRJEB22863_NSCLC_msv_dist_std.RData")
PRJEB22863_NSCLC_msv_dist_std<-msv_dist_std
load("01.cleanData/SV_all/distMat/PRJEB22863_RCC_msv_dist_std.RData")
PRJEB22863_RCC_msv_dist_std<-msv_dist_std
load("01.cleanData/SV_all/distMat/PRJNA751792_msv_dist_std.RData")
PRJNA751792_msv_dist_std<-msv_dist_std
load("01.cleanData/SV_all/distMat/PRJNA541981_msv_dist_std.RData")
PRJNA541981_msv_dist_std<-msv_dist_std


######################### 2 Species-level association
### 2.1 Prepare abundance table
if (!dir.exists("06.Species_genetic_association")) {dir.create("06.Species_genetic_association")}
if (!dir.exists("06.Species_genetic_association/RData")) {dir.create("06.Species_genetic_association/RData")}

###############################################################################
######## Prepare covariate table

# PRJEB22863
PRJEB22863_abun$Others<-1-rowSums(PRJEB22863_abun)
PRJEB22863_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB22863_abun)), transform="clr") %>%as.data.frame
PRJEB22863_abun_clr <- PRJEB22863_abun_clr[match(rownames(PRJEB22863_abun), rownames(PRJEB22863_abun_clr)),]
rownames(PRJEB22863_abun_clr) <- rownames(PRJEB22863_abun)
PRJEB22863_covar<-cbind(PRJEB22863_basic,PRJEB22863_abun_clr)

# PRJEB22863_NSCLC
PRJEB22863_NSCLC_abun$Others<-1-rowSums(PRJEB22863_NSCLC_abun)
PRJEB22863_NSCLC_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB22863_NSCLC_abun)), transform="clr") %>%as.data.frame
PRJEB22863_NSCLC_abun_clr <- PRJEB22863_NSCLC_abun_clr[match(rownames(PRJEB22863_NSCLC_abun), rownames(PRJEB22863_NSCLC_abun_clr)),]
rownames(PRJEB22863_NSCLC_abun_clr) <- rownames(PRJEB22863_NSCLC_abun)
PRJEB22863_NSCLC_covar<-cbind(PRJEB22863_NSCLC_basic,PRJEB22863_NSCLC_abun_clr)

# PRJEB22863_RCC
PRJEB22863_RCC_abun$Others<-1-rowSums(PRJEB22863_RCC_abun)
PRJEB22863_RCC_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB22863_RCC_abun)), transform="clr") %>%as.data.frame
PRJEB22863_RCC_abun_clr <- PRJEB22863_RCC_abun_clr[match(rownames(PRJEB22863_RCC_abun), rownames(PRJEB22863_RCC_abun_clr)),]
rownames(PRJEB22863_RCC_abun_clr) <- rownames(PRJEB22863_RCC_abun)
PRJEB22863_RCC_covar<-cbind(PRJEB22863_RCC_basic,PRJEB22863_RCC_abun_clr)


# PRJNA397906
PRJNA397906_abun$Others<-1-rowSums(PRJNA397906_abun)
PRJNA397906_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA397906_abun)), transform="clr") %>%as.data.frame
PRJNA397906_abun_clr <- PRJNA397906_abun_clr[match(rownames(PRJNA397906_abun), rownames(PRJNA397906_abun_clr)),]
rownames(PRJNA397906_abun_clr) <- rownames(PRJNA397906_abun)
PRJNA397906_covar<-cbind(PRJNA397906_basic,PRJNA397906_abun_clr)

# PRJNA751792
PRJNA751792_abun$Others<-1-rowSums(PRJNA751792_abun)
PRJNA751792_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA751792_abun)), transform="clr") %>%as.data.frame
PRJNA751792_abun_clr <- PRJNA751792_abun_clr[match(rownames(PRJNA751792_abun), rownames(PRJNA751792_abun_clr)),]
rownames(PRJNA751792_abun_clr) <- rownames(PRJNA751792_abun)
PRJNA751792_covar<-cbind(PRJNA751792_basic,PRJNA751792_abun_clr)

# PRJNA762360
PRJNA762360_abun$Others<-1-rowSums(PRJNA762360_abun)
PRJNA762360_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA762360_abun)), transform="clr") %>%as.data.frame
PRJNA762360_abun_clr <- PRJNA762360_abun_clr[match(rownames(PRJNA762360_abun), rownames(PRJNA762360_abun_clr)),]
rownames(PRJNA762360_abun_clr) <- rownames(PRJNA762360_abun)
PRJNA762360_covar<-cbind(PRJNA762360_basic,PRJNA762360_abun_clr)


# PRJNA770295
PRJNA770295_abun$Others<-1-rowSums(PRJNA770295_abun)
PRJNA770295_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA770295_abun)), transform="clr") %>%as.data.frame
PRJNA770295_abun_clr <- PRJNA770295_abun_clr[match(rownames(PRJNA770295_abun), rownames(PRJNA770295_abun_clr)),]
rownames(PRJNA770295_abun_clr) <- rownames(PRJNA770295_abun)
PRJNA770295_covar<-cbind(PRJNA770295_basic,PRJNA770295_abun_clr)


# PRJEB43119
PRJEB43119_abun$Others<-1-rowSums(PRJEB43119_abun)
PRJEB43119_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB43119_abun)), transform="clr") %>%as.data.frame
PRJEB43119_abun_clr <- PRJEB43119_abun_clr[match(rownames(PRJEB43119_abun), rownames(PRJEB43119_abun_clr)),]
rownames(PRJEB43119_abun_clr) <- rownames(PRJEB43119_abun)
PRJEB43119_covar<-cbind(PRJEB43119_basic,PRJEB43119_abun_clr)

# PRJNA541981
PRJNA541981_abun$Others<-1-rowSums(PRJNA541981_abun)
PRJNA541981_abun_clr<-abundances(x=as.data.frame(na.omit(PRJNA541981_abun)), transform="clr") %>%as.data.frame
PRJNA541981_abun_clr <- PRJNA541981_abun_clr[match(rownames(PRJNA541981_abun), rownames(PRJNA541981_abun_clr)),]
rownames(PRJNA541981_abun_clr) <- rownames(PRJNA541981_abun)
PRJNA541981_covar<-cbind(PRJNA541981_basic,PRJNA541981_abun_clr)

### 2.2 Association between genetic makeup and prognosis 
### covar1 是需要校正的协变量名称
### all_Prog 是因变量（预后相关的变量）
### all_covar是包含自变量的数据集 
### info是物种的信息集

covar1 <- c("age_bins","gender_num","cohort")
covar2 <- c("age_bins","gender_num")

### pfs有些数据集只有一套数据有信息，如果校正cohort，取值只有1个，不合理
### different datasets  (全部为缺失值的变量不能算进去）

covar_PRJEB22863<-PRJEB22863_covar[,c("response_code","pfs_12_months")]
PRJEB22863_adonis_res <- my_adonis_terms_noadjAbun(PRJEB22863_msv_dist_std, covar_PRJEB22863, PRJEB22863_covar,covar2, France_info)
save(PRJEB22863_adonis_res, file = "06.Species_genetic_association/RData/PRJEB22863_adonis_res.RData")

covar_PRJEB22863_NSCLC<-PRJEB22863_NSCLC_covar[,c("response_code","pfs_12_months")]
PRJEB22863_NSCLC_adonis_res <- my_adonis_terms_noadjAbun(PRJEB22863_NSCLC_msv_dist_std, covar_PRJEB22863_NSCLC, PRJEB22863_NSCLC_covar,covar2, France_info)
save(PRJEB22863_NSCLC_adonis_res, file = "06.Species_genetic_association/RData/PRJEB22863_NSCLC_adonis_res.RData")

covar_PRJEB22863_NSCLC<-PRJEB22863_NSCLC_covar[,c("response_code","pfs_12_months")]
PRJEB22863_NSCLC_adonis_res <- my_adonis_terms_noadjAbun(PRJEB22863_NSCLC_msv_dist_std, covar_PRJEB22863_NSCLC, PRJEB22863_NSCLC_covar,covar2, France_info)
save(PRJEB22863_NSCLC_adonis_res, file = "06.Species_genetic_association/RData/PRJEB22863_NSCLC_adonis_res_noadjust.RData")

covar_PRJEB22863_RCC<-PRJEB22863_RCC_covar[,c("response_code","pfs_12_months")]
PRJEB22863_RCC_adonis_res <- my_adonis_terms_noadjAbun(PRJEB22863_RCC_msv_dist_std, covar_PRJEB22863_RCC, PRJEB22863_RCC_covar,covar2, France_info)
save(PRJEB22863_RCC_adonis_res, file = "06.Species_genetic_association/RData/PRJEB22863_RCC_adonis_res.RData")

covar_PRJNA751792<-PRJNA751792_covar[,c("response_code","response_code")]
PRJNA751792_adonis_res <- my_adonis_terms_noadjAbun(PRJNA751792_msv_dist_std, covar_PRJNA751792, PRJNA751792_covar, covar2, France_info)
save(PRJNA751792_adonis_res, file = "06.Species_genetic_association/RData/PRJNA751792_adonis_res.RData")

covar_PRJNA397906<-PRJNA397906_covar[,c("response_code","response_code")]
PRJNA397906_adonis_res <- my_adonis_terms_noadjAbun(PRJNA397906_msv_dist_std, covar_PRJNA397906, PRJNA397906_covar, covar2, USA_UK_info)
save(PRJNA397906_adonis_res, file = "06.Species_genetic_association/RData/PRJNA397906_adonis_res.RData")

covar_PRJNA762360<-PRJNA762360_covar[,c("response_code","irAEs","pfs_12_months")]
PRJNA762360_adonis_res <- my_adonis_terms_noadjAbun(PRJNA762360_msv_dist_std, covar_PRJNA762360, PRJNA762360_covar, covar2, USA_UK_info)
save(PRJNA762360_adonis_res, file = "06.Species_genetic_association/RData/PRJNA762360_adonis_res.RData")

covar_PRJNA770295<-PRJNA770295_covar[,c("response_code","pfs_12_months")]
PRJNA770295_adonis_res <- my_adonis_terms_noadjAbun(PRJNA770295_msv_dist_std, covar_PRJNA770295, PRJNA770295_covar, covar2, USA_UK_info)
save(PRJNA770295_adonis_res, file = "06.Species_genetic_association/RData/PRJNA770295_adonis_res.RData")

covar_PRJEB43119<-PRJEB43119_covar[,c("response_code","pfs_12_months")]
PRJEB43119_adonis_res <- my_adonis_terms_noadjAbun(PRJEB43119_msv_dist_std, covar_PRJEB43119, PRJEB43119_covar, covar2, USA_UK_info)
save(PRJEB43119_adonis_res, file = "06.Species_genetic_association/RData/PRJEB43119_adonis_res.RData")

covar_PRJNA541981<-PRJNA541981_covar[,c("pfs_12_months","pfs_12_months")]
PRJNA541981_adonis_res <- my_adonis_terms_noadjAbun(PRJNA541981_msv_dist_std, covar_PRJNA541981, PRJNA541981_covar, covar2, USA_UK_info)
save(PRJNA541981_adonis_res, file = "06.Species_genetic_association/RData/PRJNA541981_adonis_res.RData")


## meta-analysis
#load("06.Species_genetic_association/RData/all_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJEB22863_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJEB22863_NSCLC_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJEB22863_RCC_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJNA397906_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJNA751792_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJNA762360_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJNA770295_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJEB43119_adonis_res.RData")
load("06.Species_genetic_association/RData/PRJNA541981_adonis_res.RData")


#all_adonis_res.table <- all_adonis_res$table
PRJEB22863_adonis_res.table <- PRJEB22863_adonis_res$table
PRJEB22863_NSCLC_adonis_res.table <- PRJEB22863_NSCLC_adonis_res$table
PRJEB22863_RCC_adonis_res.table <- PRJEB22863_RCC_adonis_res$table
PRJNA397906_adonis_res.table <- PRJNA397906_adonis_res$table
PRJNA751792_adonis_res.table  <- PRJNA751792_adonis_res$table
PRJNA762360_adonis_res.table  <- PRJNA762360_adonis_res$table
PRJNA770295_adonis_res.table  <- PRJNA770295_adonis_res$table
PRJEB43119_adonis_res.table  <- PRJEB43119_adonis_res$table
PRJNA541981_adonis_res.table  <- PRJNA541981_adonis_res$table

PRJEB22863_adonis_resp<-subset(PRJEB22863_adonis_res.table,Prog=="response_code",drop=T)
PRJEB22863_adonis_resp<-PRJEB22863_adonis_resp[,-c(1,2)]
colnames(PRJEB22863_adonis_res.table)<-c(paste("RoutyB_2018.",colnames(PRJEB22863_adonis_resp)[1:5],sep = ""))

PRJEB22863_NSCLC_adonis_resp<-subset(PRJEB22863_NSCLC_adonis_res.table,Prog=="response_code",drop=T)
PRJEB22863_NSCLC_adonis_resp<-PRJEB22863_NSCLC_adonis_resp[,-c(1,2)]
colnames(PRJEB22863_NSCLC_adonis_resp)<-c(paste("RoutyB_2018_NSCLC.",colnames(PRJEB22863_NSCLC_adonis_resp)[1:5],sep = ""))

PRJEB22863_RCC_adonis_resp<-subset(PRJEB22863_RCC_adonis_res.table,Prog=="response_code",drop=T)
colnames(PRJEB22863_RCC_adonis_resp)<-c(paste("RoutyB_2018_RCC.",colnames(PRJEB22863_RCC_adonis_resp)[1:7],sep = ""))

PRJNA397906_adonis_resp<-subset(PRJNA397906_adonis_res.table,Prog=="response_code",drop=T)
PRJNA397906_adonis_resp<-PRJNA397906_adonis_resp[,-c(1,2)]
colnames(PRJNA397906_adonis_resp)<-c(paste("FrankelAE_2017.",colnames(PRJNA397906_adonis_resp)[1:5],sep = ""))

PRJNA751792_adonis_resp<-subset(PRJNA751792_adonis_res.table,Prog=="response_code",drop=T)
##PRJNA751792_adonis_resp<-PRJNA751792_adonis_resp[,-c(1,2)]
colnames(PRJNA751792_adonis_resp)<-c(paste("DerosaL_2022.",colnames(PRJNA751792_adonis_resp)[1:7],sep = ""))

PRJNA762360_adonis_resp<-subset(PRJNA762360_adonis_res.table,Prog=="response_code",drop=T)
##PRJNA762360_adonis_resp<-PRJNA762360_adonis_resp[,-c(1,2)]
colnames(PRJNA762360_adonis_resp)<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_adonis_resp)[1:7],sep = ""))

PRJNA770295_adonis_resp<-subset(PRJNA770295_adonis_res.table,Prog=="response_code",drop=T)
PRJNA770295_adonis_resp<-PRJNA770295_adonis_resp[,-c(1,2)]
colnames(PRJNA770295_adonis_resp)<-c(paste("SpencerCN_2021.",colnames(PRJNA770295_adonis_resp)[1:5],sep = ""))

PRJEB43119_adonis_resp<-subset(PRJEB43119_adonis_res.table,Prog=="response_code",drop=T)
PRJEB43119_adonis_resp<-PRJEB43119_adonis_resp[,-c(1,2)]
colnames(PRJEB43119_adonis_resp)<-c(paste("LeeKA_2022.",colnames(PRJEB43119_adonis_resp)[1:5],sep = ""))


#################################################################################################
### response 
### mel 
mel_cbind_adonis_resp<-cbind(PRJNA762360_adonis_resp,PRJNA397906_adonis_resp,
PRJNA770295_adonis_resp,PRJEB43119_adonis_resp)
mel_resp_adonis <- my_batch_meta_p(mel_cbind_adonis_resp, c("FrankelAE_2017","McCullochJA_2022",
"SpencerCN_2021","LeeKA_2022"),seq(5,22,by=5),seq(4,22,by=5))$table
write.csv(mel_resp_adonis, "06.Species_genetic_association/mel_resp_adonis.csv",row.names=F)

### nsclc 
nsclc_cbind_adonis_resp<-cbind(PRJNA751792_adonis_resp,PRJEB22863_NSCLC_adonis_resp)
nsclc_resp_adonis <- my_batch_meta_p(nsclc_cbind_adonis_resp, c("DerosaL_2022","RoutyB_2018_NSCLC"),
seq(5,12,by=5),seq(4,12,by=5))$table
write.csv(nsclc_resp_adonis, "06.Species_genetic_association/nsclc_resp_adonis.csv",row.names=F)

### rcc
PRJEB22863_RCC_adonis_resp<-subset(PRJEB22863_RCC_adonis_res.table,Prog=="response_code",drop=T)
colnames(PRJEB22863_RCC_adonis_resp)<-c(paste("RoutyB_2018_RCC.",colnames(PRJEB22863_RCC_adonis_resp)[1:7],sep = ""))
write.csv(PRJEB22863_RCC_adonis_resp, "06.Species_genetic_association/rcc_resp_adonis.csv",row.names=F)


#####################################################################
##  irAEs
PRJNA762360_adonis_irAEs<-subset(PRJNA762360_adonis_res.table,Prog=="irAEs",drop=T)
colnames(PRJNA762360_adonis_irAEs)<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_adonis_irAEs)[1:7],sep = ""))
write.csv(PRJNA762360_adonis_irAEs, "06.Species_genetic_association/mel_irAEs_adonis.csv",row.names=F)


##########################################################################################################
### pfs  12 months and species

PRJNA762360_adonis_pfs_12<-subset(PRJNA762360_adonis_res.table,Prog=="pfs_12_months",drop=T)
colnames(PRJNA762360_adonis_pfs_12)<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_adonis_pfs_12)[1:7],sep = ""))

PRJNA770295_adonis_pfs_12<-subset(PRJNA770295_adonis_res.table,Prog=="pfs_12_months",drop=T)[,-c(1,2)]
colnames(PRJNA770295_adonis_pfs_12)<-c(paste("SpencerCN_2021.",colnames(PRJNA770295_adonis_pfs_12)[1:5],sep = ""))

PRJEB43119_adonis_pfs_12<-subset(PRJEB43119_adonis_res.table,Prog=="pfs_12_months",drop=T)[,-c(1,2)]
colnames(PRJEB43119_adonis_pfs_12)<-c(paste("LeeKA_2022.",colnames(PRJEB43119_adonis_pfs_12)[1:5],sep = ""))

PRJNA541981_adonis_pfs_12<-subset(PRJNA541981_adonis_res.table,Prog=="pfs_12_months",drop=T)[,-c(1,2)]
colnames(PRJNA541981_adonis_pfs_12)<-c(paste("PetersBA_2019.",colnames(PRJNA541981_adonis_pfs_12)[1:5],sep = ""))

mel_cbind_adonis_pfs_12<-cbind(PRJNA762360_adonis_pfs_12,PRJNA770295_adonis_pfs_12,PRJEB43119_adonis_pfs_12,PRJNA541981_adonis_pfs_12)

mel_adonis_res <- my_batch_meta_p(mel_cbind_adonis_pfs_12, c("McCullochJA_2022","SpencerCN_2021","LeeKA_2022","PetersBA_2019"),
seq(5,22,by=5),seq(4,22,by=5))$table

write.csv(mel_adonis_res, "06.Species_genetic_association/mel_pfs_12_months_adonis.csv",row.names=F)


#############################################################
###  2.3 Association between abundance and prognosis
### overall 

covar1 <- c("age_bins","gender_num","cohort")
covar2 <- c("age_bins","gender_num")
nc_abun_France<-ncol(PRJEB22863_abun_clr)
nc_abun_USA_UK<-ncol(PRJNA770295_abun_clr)

covar_PRJEB22863<-PRJEB22863_covar[,c("response_code","pfs_12_months")]
PRJEB22863_sv_abun_lm_res <- lm_btw_mats(covar_PRJEB22863,PRJEB22863_abun_clr[,-nc_abun_France],PRJEB22863_basic,covar2)
save(PRJEB22863_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJEB22863_abun_res.RData")

covar_PRJEB22863_NSCLC<-PRJEB22863_NSCLC_covar[,c("response_code","pfs_12_months")]
PRJEB22863_NSCLC_sv_abun_lm_res <- lm_btw_mats(covar_PRJEB22863_NSCLC,PRJEB22863_NSCLC_abun_clr[,-nc_abun_France],PRJEB22863_NSCLC_basic,covar2)
save(PRJEB22863_NSCLC_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJEB22863_NSCLC_abun_res.RData")

covar_PRJEB22863_RCC<-PRJEB22863_RCC_covar[,c("response_code","pfs_12_months")]
PRJEB22863_RCC_sv_abun_lm_res <- lm_btw_mats(covar_PRJEB22863_RCC,PRJEB22863_RCC_abun_clr[,-nc_abun_France],PRJEB22863_RCC_basic,covar2)
save(PRJEB22863_RCC_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJEB22863_RCC_abun_res.RData")

covar_PRJNA751792<-PRJNA751792_covar[,c("response_code","response_code")]
PRJNA751792_sv_abun_lm_res <- lm_btw_mats(covar_PRJNA751792,PRJNA751792_abun_clr[,-nc_abun_France],PRJNA751792_basic,covar2)
save(PRJNA751792_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJNA751792_abun_res.RData")

covar_PRJNA397906<-PRJNA397906_covar[,c("response_code","response_code")]
PRJNA397906_sv_abun_lm_res <- lm_btw_mats(covar_PRJNA397906,PRJNA397906_abun_clr[,-nc_abun_USA_UK],PRJNA397906_basic,covar2)
save(PRJNA397906_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJNA397906_abun_res.RData")

covar_PRJNA762360<-PRJNA762360_covar[,c("response_code","irAEs","pfs_12_months")]
PRJNA762360_sv_abun_lm_res <- lm_btw_mats(covar_PRJNA762360,PRJNA762360_abun_clr[,-nc_abun_USA_UK],PRJNA762360_basic,covar2)
save(PRJNA762360_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJNA762360_abun_res.RData")

covar_PRJNA770295<-PRJNA770295_covar[,c("response_code","pfs_12_months")]
PRJNA770295_sv_abun_lm_res <- lm_btw_mats(covar_PRJNA770295,PRJNA770295_abun_clr[,-nc_abun_USA_UK],PRJNA770295_basic,covar2)
save(PRJNA770295_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJNA770295_abun_res.RData")

covar_PRJEB43119<-PRJEB43119_covar[,c("response_code","pfs_12_months")]
PRJEB43119_sv_abun_lm_res <- lm_btw_mats(covar_PRJEB43119,PRJEB43119_abun_clr[,-nc_abun_USA_UK],PRJEB43119_basic,covar2)
save(PRJEB43119_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJEB43119_abun_res.RData")

covar_PRJNA541981<-PRJNA541981_covar[,c("pfs_12_months","pfs_12_months")]
PRJNA541981_sv_abun_lm_res <- lm_btw_mats(covar_PRJNA541981,PRJNA541981_abun_clr[,-nc_abun_USA_UK],PRJNA541981_basic,covar2)
save(PRJNA541981_sv_abun_lm_res, file = "06.Species_genetic_association/RData/PRJNA541981_abun_res.RData")


###########################################################################################################
### meta-analysis

load("06.Species_genetic_association/RData/PRJEB22863_abun_res.RData")
load("06.Species_genetic_association/RData/PRJEB22863_NSCLC_abun_res.RData")
load("06.Species_genetic_association/RData/PRJEB22863_RCC_abun_res.RData")
load("06.Species_genetic_association/RData/PRJNA397906_abun_res.RData")
load("06.Species_genetic_association/RData/PRJNA751792_abun_res.RData")
load("06.Species_genetic_association/RData/PRJNA762360_abun_res.RData")
load("06.Species_genetic_association/RData/PRJNA770295_abun_res.RData")
load("06.Species_genetic_association/RData/PRJEB43119_abun_res.RData")
load("06.Species_genetic_association/RData/PRJNA541981_abun_res.RData")

PRJEB22863_sv_abun.table <- PRJEB22863_sv_abun_lm_res
PRJEB22863_NSCLC_sv_abun.table <- PRJEB22863_NSCLC_sv_abun_lm_res
PRJEB22863_RCC_sv_abun.table <- PRJEB22863_RCC_sv_abun_lm_res
PRJNA397906_sv_abun.table <- PRJNA397906_sv_abun_lm_res
PRJNA751792_sv_abun.table  <- PRJNA751792_sv_abun_lm_res
PRJNA762360_sv_abun.table  <- PRJNA762360_sv_abun_lm_res
PRJNA770295_sv_abun.table  <- PRJNA770295_sv_abun_lm_res
PRJEB43119_sv_abun.table  <- PRJEB43119_sv_abun_lm_res
PRJNA541981_sv_abun.table  <- PRJNA541981_sv_abun_lm_res

PRJEB22863_sv_abun_resp<-subset(PRJEB22863_sv_abun.table,Y=="response_code",drop=T)
PRJEB22863_sv_abun_resp<-PRJEB22863_sv_abun_resp[,-c(1,2)]
colnames(PRJEB22863_sv_abun_resp)<-c(paste("RoutyB_2018.",colnames(PRJEB22863_sv_abun_resp)[1:6],sep = ""))

PRJEB22863_NSCLC_sv_abun_resp<-subset(PRJEB22863_sv_abun.table,Y=="response_code",drop=T)
PRJEB22863_NSCLC_sv_abun_resp<-PRJEB22863_NSCLC_sv_abun_resp[,-c(1,2)]
colnames(PRJEB22863_NSCLC_sv_abun_resp)<-c(paste("RoutyB_2018_NSCLC.",colnames(PRJEB22863_NSCLC_sv_abun_resp)[1:6],sep = ""))

PRJEB22863_RCC_sv_abun_resp<-subset(PRJEB22863_sv_abun.table,Y=="response_code",drop=T)
colnames(PRJEB22863_RCC_sv_abun_resp)<-c(paste("RoutyB_2018_RCC.",colnames(PRJEB22863_RCC_sv_abun_resp)[1:8],sep = ""))

PRJNA397906_sv_abun_resp<-subset(PRJNA397906_sv_abun.table,Y=="response_code",drop=T)
PRJNA397906_sv_abun_resp<-PRJNA397906_sv_abun_resp[,-c(1,2)]
colnames(PRJNA397906_sv_abun_resp)<-c(paste("FrankelAE_2017.",colnames(PRJNA397906_sv_abun_resp)[1:6],sep = ""))

PRJNA751792_sv_abun_resp<-subset(PRJNA751792_sv_abun.table,Y=="response_code",drop=T)
##PRJNA751792_sv_abun_resp<-PRJNA751792_sv_abun_resp[,-c(1,2)]
colnames(PRJNA751792_sv_abun_resp)<-c(paste("DerosaL_2022.",colnames(PRJNA751792_sv_abun_resp)[1:8],sep = ""))

PRJNA762360_sv_abun_resp<-subset(PRJNA762360_sv_abun.table,Y=="response_code",drop=T)
##PRJNA762360_sv_abun_resp<-PRJNA762360_sv_abun_resp[,-c(1,2)]
colnames(PRJNA762360_sv_abun_resp)<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_sv_abun_resp)[1:8],sep = ""))

PRJNA770295_sv_abun_resp<-subset(PRJNA770295_sv_abun.table,Y=="response_code",drop=T)
PRJNA770295_sv_abun_resp<-PRJNA770295_sv_abun_resp[,-c(1,2)]
colnames(PRJNA770295_sv_abun_resp)<-c(paste("SpencerCN_2021.",colnames(PRJNA770295_sv_abun_resp)[1:6],sep = ""))

PRJEB43119_sv_abun_resp<-subset(PRJEB43119_sv_abun.table,Y=="response_code",drop=T)
PRJEB43119_sv_abun_resp<-PRJEB43119_sv_abun_resp[,-c(1,2)]
colnames(PRJEB43119_sv_abun_resp)<-c(paste("LeeKA_2022.",colnames(PRJEB43119_sv_abun_resp)[1:6],sep = ""))


#################################################################################################
### response 
### mel 
mel_cbind_sv_abun_resp<-cbind(PRJNA762360_sv_abun_resp,PRJNA397906_sv_abun_resp,
PRJNA770295_sv_abun_resp,PRJEB43119_sv_abun_resp)
my_batch_meta_p_lr(mel_cbind_sv_abun_resp,seq(3,26,6),seq(4,26,6),"06.Species_genetic_association/mel_resp_abun.csv")

### nsclc 
nsclc_cbind_sv_abun_resp<-cbind(PRJNA751792_sv_abun_resp,PRJEB22863_NSCLC_sv_abun_resp)
my_batch_meta_p_lr(nsclc_cbind_sv_abun_resp,seq(3,14,6),seq(4,14,6),"06.Species_genetic_association/nsclc_resp_abun.csv")

#### rcc 
PRJEB22863_RCC_abun_resp<-subset(PRJEB22863_RCC_sv_abun.table,Y=="response_code",drop=T)
colnames(PRJEB22863_RCC_abun_resp)<-c(paste("RoutyB_2018_RCC.",colnames(PRJEB22863_RCC_abun_resp)[1:8],sep = ""))
write.csv(PRJEB22863_RCC_abun_resp, "06.Species_genetic_association/rcc_resp_abun.csv",row.names=F)

##########################################################################################################
### pfs  12 months and species

PRJNA762360_abun_pfs_12<-subset(PRJNA762360_sv_abun.table,Y=="pfs_12_months",drop=T)
colnames(PRJNA762360_abun_pfs_12)<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_abun_pfs_12)[1:8],sep = ""))

PRJNA770295_abun_pfs_12<-subset(PRJNA770295_sv_abun.table,Y=="pfs_12_months",drop=T)[,-c(1,2)]
colnames(PRJNA770295_abun_pfs_12)<-c(paste("SpencerCN_2021.",colnames(PRJNA770295_abun_pfs_12)[1:6],sep = ""))

PRJEB43119_abun_pfs_12<-subset(PRJEB43119_sv_abun.table,Y=="pfs_12_months",drop=T)[,-c(1,2)]
colnames(PRJEB43119_abun_pfs_12)<-c(paste("LeeKA_2022.",colnames(PRJEB43119_abun_pfs_12)[1:6],sep = ""))

PRJNA541981_abun_pfs_12<-subset(PRJNA541981_sv_abun.table,Y=="pfs_12_months",drop=T)[,-c(1,2)]
colnames(PRJNA541981_abun_pfs_12)<-c(paste("PetersBA_2019.",colnames(PRJNA541981_abun_pfs_12)[1:6],sep = ""))

mel_cbind_abun_pfs_12<-cbind(PRJNA762360_abun_pfs_12,PRJNA770295_abun_pfs_12,PRJEB43119_abun_pfs_12,PRJNA541981_abun_pfs_12)
my_batch_meta_p_lr(mel_cbind_abun_pfs_12,seq(3,26,6),seq(4,26,6),"06.Species_genetic_association/mel_pfs_12_months_abun.csv")

#### irAEs
PRJNA762360_abun_irAEs<-subset(PRJNA762360_sv_abun.table,Y=="irAEs",drop=T)
colnames(PRJNA762360_abun_irAEs)<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_abun_irAEs)[1:8],sep = ""))
write.csv(PRJNA762360_abun_irAEs, "06.Species_genetic_association/mel_irAEs_abun.csv",row.names=F)


######################################################################################################
### 3 Visualization   (Mel)
### 3.1 Preparation
### SVs association

mel_resp_adonis<-read.csv("06.Species_genetic_association/mel_resp_adonis.csv",header=T)
mel_resp_adonis<-mel_resp_adonis[,c(1,(ncol(mel_resp_adonis)-4))]
mel_resp_adonis$Prog<-"Melanoma_response"
colnames(mel_resp_adonis)<-c("Species","p","Prog")

mel_pfs_12_months_adonis<-read.csv("06.Species_genetic_association/mel_pfs_12_months_adonis.csv",header=T)
mel_pfs_12_months_adonis<-mel_pfs_12_months_adonis[,c(1,(ncol(mel_pfs_12_months_adonis)-4))]
mel_pfs_12_months_adonis$Prog<-"Melanoma_pfs_12_months"
colnames(mel_pfs_12_months_adonis)<-c("Species","p","Prog")

mel_irAEs_adonis<-read.csv("06.Species_genetic_association/mel_irAEs_adonis.csv",header=T)
mel_irAEs_adonis<-mel_irAEs_adonis[,c(1,(ncol(mel_irAEs_adonis)-3))]
mel_irAEs_adonis$Prog<-"Melanoma_irAEs"
colnames(mel_irAEs_adonis)<-c("Species","p","Prog")

all_adonis.table<-rbind(mel_resp_adonis,mel_pfs_12_months_adonis,mel_irAEs_adonis)
save(all_adonis.table, file = "06.Species_genetic_association/RData/Mel_all_adonis.table.RData")


#### abdundance significant
mel_resp_abun<-read.csv("06.Species_genetic_association/mel_resp_abun.csv",header=T)
mel_resp_abun<-mel_resp_abun[,c(2,(ncol(mel_resp_abun)-5))]
mel_resp_abun$Prog<-"Melanoma_response"
colnames(mel_resp_abun)<-c("Species","p","Prog")

mel_pfs_12_months_abun<-read.csv("06.Species_genetic_association/mel_pfs_12_months_abun.csv",header=T)
mel_pfs_12_months_abun<-mel_pfs_12_months_abun[,c(2,(ncol(mel_pfs_12_months_abun)-5))]
mel_pfs_12_months_abun$Prog<-"Melanoma_pfs_12_months"
colnames(mel_pfs_12_months_abun)<-c("Species","p","Prog")

mel_irAEs_abun<-read.csv("06.Species_genetic_association/mel_irAEs_abun.csv",header=T)
mel_irAEs_abun<-mel_irAEs_abun[,c(2,(ncol(mel_irAEs_abun)-3))]
mel_irAEs_abun$Prog<-"Melanoma_irAEs"
colnames(mel_irAEs_abun)<-c("Species","p","Prog")

all_abun.table<-rbind(mel_resp_abun,mel_pfs_12_months_abun,mel_irAEs_abun)
save(all_abun.table, file = "06.Species_genetic_association/RData/Mel_all_abun.table.RData")

all_abun.table$Species<-USA_UK_info$Short_name[match(all_abun.table$Species, USA_UK_info$organism)]

sv_assoc_id<-paste(all_adonis.table$Species, all_adonis.table$Prog,sep = "_")
abun_assoc_id<-paste(all_abun.table$Species, all_abun.table$Prog, sep = "_")

all_abun.table<-all_abun.table[match(sv_assoc_id,abun_assoc_id),]

colnames(all_abun.table)<-paste("Abun",colnames(all_abun.table),sep = ".")
save(all_abun.table, file = "06.Species_genetic_association/RData/Mel_all_abun.table.RData")

species.table<-cbind(all_adonis.table,all_abun.table)

species.table$sv.MetaSigAssoc<-rep('No', nrow(species.table))
species.table$sv.MetaSigAssoc[species.table$p< 0.05]<-'Yes'
#species.table$sv.MetaSigAssoc[species.table$fdr< 0.05]<-'Yes'
species.table$sv.MetaSigAssoc[is.na(species.table$p)==T]<-'Unknown'

species.table$abun.MetaSigAssoc<-rep('No', nrow(species.table))
species.table$abun.MetaSigAssoc[species.table$Abun.p< 0.05]<-'Yes'
#species.table$abun.MetaSigAssoc[species.table$Abun.fdr< 0.05]<-'Yes'

species.sig.table<-species.table[,] %>%
  .[.$sv.MetaSigAssoc=="Yes" | .$abun.MetaSigAssoc == "Yes",]

write.table(species.sig.table, "06.Species_genetic_association/Mel_species.sig.tsv",sep = "\t", 
   col.names = T, row.names = F, quote = F)

write.csv(species.sig.table, "06.Species_genetic_association/Mel_species.sig.csv",row.names = F)

unique(species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]$Species)
unique(species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]$Prog)

unique(species.sig.table[species.sig.table$abun.MetaSigAssoc=="Yes",]$Species)
unique(species.sig.table[species.sig.table$abun.MetaSigAssoc=="Yes",]$Prog)

species.sv.sig.table<-species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]

#### more than two datasets
#cor.test(species_Prog.sv.sig.table$PRJEB22863.R2,species_Prog.sv.sig.table$PRJNA397906.R2)
#plot(species_Prog.sv.sig.table$PRJEB22863.R2,species_Prog.sv.sig.table$PRJNA397906.R2)

### 3.4 Venn diagram
species_count<-table(species.table$sv.MetaSigAssoc,species.table$abun.MetaSigAssoc)

#pdf("06.Species_genetic_association/species_Prog_count.venn.pdf", width = 2, height = 2)
tiff(file = "pics/FS2C_Mel_species_count.venn.tiff", width =600, height =600, res =300) 

draw.pairwise.venn(species_count[3,1]+species_count[3,2],
                   species_count[1,2]+species_count[2,2]+species_count[3,2],
                   species_count[3,2], 
                   category = c("Genetics", "Abundance"), lty = rep("blank",2), 
                   fill =c("#4472c4", "#00b050"), alpha = rep(0.5, 2), cat.pos = c(0, 0), 
                   cat.cex = c(0.5, 0.5),cat.dist = rep(0.025, 2), scaled = F)

dev.off()
 

##############################################################################
###  3.5 heatmap
load("06.Species_genetic_association/RData/Mel_all_adonis.table.Rdata")
load("06.Species_genetic_association/RData/Mel_all_abun.table.RData")

## matrix
adonis_res.P<-all_adonis.table[,c("Prog", "Species", "p")] %>%
  spread("Prog", "p")
adonis_res.P<-data.frame(adonis_res.P, row.names = "Species")

all_abun_res.P<-all_abun.table[,c("Abun.Prog", "Abun.Species", "Abun.p")] %>%
  spread("Abun.Prog", "Abun.p")

all_abun_res.P<-data.frame(all_abun_res.P, row.names = "Abun.Species")
all_abun_res.P<-all_abun_res.P[match(rownames(adonis_res.P),rownames(all_abun_res.P)),
                                                     match(colnames(adonis_res.P),colnames(all_abun_res.P))]

# color
species.color<-matrix(0, nrow = nrow(adonis_res.P), ncol = ncol(adonis_res.P))
rownames(species.color)<-rownames(adonis_res.P)
colnames(species.color)<-colnames(adonis_res.P)

species.color[adonis_res.P<0.05&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 1        #"purple", both assoc
species.color[adonis_res.P>=0.05&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 2       #"yellow", abundance assoc
species.color[adonis_res.P<0.05&all_abun_res.P>=0.05&!is.na(all_abun_res.P)]<- 3       #"blue", genetic assoc
species.color[is.na(adonis_res.P)==T&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 4   #"grey", none assoc


### association number 
species.asso_number<-matrix(0, nrow = nrow(adonis_res.P), ncol = ncol(adonis_res.P))
rownames(species.asso_number)<-rownames(adonis_res.P)
colnames(species.asso_number)<-colnames(adonis_res.P)

species.asso_number[adonis_res.P<0.05&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 2        #both assoc
species.asso_number[adonis_res.P>=0.05&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 1       #abundance assoc
species.asso_number[adonis_res.P<0.05&all_abun_res.P>=0.05&!is.na(all_abun_res.P)]<- 1       #genetic assoc
species.asso_number[is.na(adonis_res.P)==T&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 0   #none assoc

species.asso<-apply(species.asso_number,1,sum)
species.asso<-species.asso[which(species.asso>0)]
species.asso<-species.asso[order(species.asso,na.last = TRUE, decreasing = TRUE)]
species_order<-names(species.asso)

## get plot tables
##species.sig.table<-species.table[species.table$fdr<0.05 | species.table$Abun.fdr<0.05,]
species.sig.table<-species.table[species.table$p<0.05 | species.table$Abun.p<0.05,]
species.sig.table<-subset(species.sig.table,is.na(sv.MetaSigAssoc)==F,drop=T)

##write.csv(species.color,"test.csv")
species.plot<- species.sig.table$Prog %>%
  as.character(.) %>%
  na.omit %>%
  .[!duplicated(.)]

species.plot.spe <- species.sig.table$Species %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

species.color.plot    <- species.color %>%
  .[match(species.plot.spe,rownames(.)),match(species.plot,colnames(.))]

colnames(species.color.plot)<-c("Response(Melanoma)","pfs >= 12 months(Melanoma)","irAEs(Melanoma)")
species.color.plot<-species.color.plot[species_order,]

tiff(file = "pics/F4A_Mel_specie_prognosis.heatmap_N50.tiff", width =1200, height =2800, res =300) 

heatmap.2(species.color.plot, 
          col=colorRampPalette(c("#e9e9e9","#993399","#E18727FF", "#0072B5FF","#666666"))(5), # white,red,yellow,blue
          #col=colorRampPalette(c("#e9e9e9","#993399","#E18727FF", "#0072B5FF"))(4),
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",
          cexCol =1.5, srtCol = 58, cexRow =1.5,
          labRow=as.expression(lapply(rownames(species.color.plot), function(a) bquote(italic(.(a))))),
          colsep = c(1:(ncol(species.color.plot)-1)), rowsep = c(1:(nrow(species.color.plot)-1)),
          sepcolor="#30303a",sepwidth = c(0.02,0.02),
          key = F,
          lmat=rbind( c(4, 4, 3), c(2, 1, 0 )), lhei = c(0.1, 3.5),lwid=c(3, 7, 4.2 ),key.title = NA,
          margins=c(16,5))

#legend('bottomright', c("Both assocation","Abundance association only","SV association only","No association","not available"), 
#box.lty =0,pch=15,
#col=c("#993399","#E18727FF", "#0072B5FF","#e9e9e9","#666666"), cex=1)

dev.off()


#######################################################################################################
### 3.6 Species pcoa  (特定的species和预后表型之间关系的展示）
## PFS 
clin_PRJNA770295<-subset(USA_UK_basic,dataset=="PRJNA770295",drop=T)
species_short_name<-"O.splanchnicus"
species_dis<-PRJNA770295_msv_dist_std[[paste("msv_",USA_UK_info$organism[match(species_short_name,USA_UK_info$Short_name)],sep = "")]]

prog_name<-"pfs_12_months"
prog_vec<-clin_PRJNA770295[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"PFS < 12 months"
prog_vec_input[prog_vec_input==1]<-"PFS >= 12 months"

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ggtitle("Melanoma (SpencerCN_2021):\nO.splanchnicus\nadonis R2=0.11;P-value:0.003")+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#E18727FF","#20854EFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 10,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

tiff(file = "pics/F4C_O.splanchnicus_pfs.tiff", width =1800, height =2300, res =600) 
print(p_spe_pcoa)
dev.off()


######################################################################
## PFS 
clin_PRJEB43119<-subset(USA_UK_basic,dataset=="PRJEB43119",drop=T)
species_short_name<-"R.gnavus"
species_dis<-PRJEB43119_msv_dist_std[[paste("msv_",USA_UK_info$organism[match(species_short_name,USA_UK_info$Short_name)],sep = "")]]

prog_name<-"pfs_12_months"
prog_vec<-clin_PRJEB43119[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"PFS < 12 months"
prog_vec_input[prog_vec_input==1]<-"PFS >= 12 months"

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ggtitle("Melanoma (LeeKA_2022):\nR.gnavus\nadonis R2=0.05;P-value:0.054")+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#E18727FF","#20854EFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 10,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

tiff(file = "pics/F4C_R.gnavus_pfs.tiff", width =1800, height =2300, res =600) 
print(p_spe_pcoa)
dev.off()


#######################################################################################
###  response PRJEB43119

species_short_name<-"D.formicigenerans"
species_dis<-PRJEB43119_msv_dist_std[[paste("msv_",USA_UK_info$organism[match(species_short_name,USA_UK_info$Short_name)],sep = "")]]

prog_name<-"response_code"
clin_PRJEB43119<-subset(USA_UK_basic,dataset=="PRJEB43119",drop=T)
#prog_vec<-clin_PRJEB43119[prog_name] %>% na.omit
#prog_vec<-mel_basic[prog_name] %>% na.omit
prog_vec<-clin_PRJEB43119[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Non-responder (SD/PD)"
prog_vec_input[prog_vec_input==1]<-"Responder (CR/PR)"

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ggtitle("Melanoma (LeeKA_2022):\nD.formicigenerans\nadonis R2=0.02;P-value:0.046")+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color =prog_vec_input ) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend =F)+
  scale_color_manual(values=c("#BC3C29FF","#6F99ADFF"))+
  guides(col = guide_legend(nrow =2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 10,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))


#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F4C_D.formicigenerans_Mel_response.tiff", width =1800, height =2300, res =600) 
print(p_spe_pcoa)
dev.off()




#######################################################################################
###  response PRJEB43119

species_short_name<-"A.muciniphila"
species_dis<-PRJEB43119_msv_dist_std[[paste("msv_",USA_UK_info$organism[match(species_short_name,USA_UK_info$Short_name)],sep = "")]]

prog_name<-"response_code"
clin_PRJEB43119<-subset(all_basic,dataset=="PRJEB43119",drop=T)
#prog_vec<-clin_PRJEB43119[prog_name] %>% na.omit
#prog_vec<-mel_basic[prog_name] %>% na.omit
prog_vec<-clin_PRJEB43119[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Non-responser"
prog_vec_input[prog_vec_input==1]<-"Responser"

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ggtitle("Melanoma (LeeKA_2022):\nA.muciniphila\nadonis R2=0.03;P-value:0.048")+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#BC3C29FF","#6F99ADFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 10,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))


#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F4C_E.rectale_Mel_response.tiff", width =1800, height =2300, res =600) 
print(p_spe_pcoa)
dev.off()




######################################################################################################
### 3 Visualization  (NSCLC and RCC)
### 3.1 Preparation
### SVs association

nsclc_resp_adonis<-read.csv("06.Species_genetic_association/nsclc_resp_adonis.csv",header=T)
nsclc_resp_adonis<-nsclc_resp_adonis[,c(1,(ncol(nsclc_resp_adonis)-4))]
nsclc_resp_adonis$Prog<-"NSCLC_response"
colnames(nsclc_resp_adonis)<-c("Species","p","Prog")

rcc_resp_adonis<-read.csv("06.Species_genetic_association/rcc_resp_adonis.csv",header=T)
rcc_resp_adonis<-rcc_resp_adonis[,c(1,(ncol(rcc_resp_adonis)-3))]
rcc_resp_adonis$Prog<-"RCC_response"
colnames(rcc_resp_adonis)<-c("Species","p","Prog")

all_adonis.table<-rbind(nsclc_resp_adonis,rcc_resp_adonis)
save(all_adonis.table, file = "06.Species_genetic_association/RData/NSCLC_RCC_all_adonis.table.RData")


#### abdundance significant
nsclc_resp_abun<-read.csv("06.Species_genetic_association/nsclc_resp_abun.csv",header=T)
nsclc_resp_abun<-nsclc_resp_abun[,c(2,(ncol(nsclc_resp_abun)-5))]
nsclc_resp_abun$Prog<-"NSCLC_response"
colnames(nsclc_resp_abun)<-c("Species","p","Prog")

rcc_resp_abun<-read.csv("06.Species_genetic_association/rcc_resp_abun.csv",header=T)
rcc_resp_abun<-rcc_resp_abun[,c(2,(ncol(rcc_resp_abun)-3))]
rcc_resp_abun$Prog<-"RCC_response"
colnames(rcc_resp_abun)<-c("Species","p","Prog")

all_abun.table<-rbind(nsclc_resp_abun,rcc_resp_abun)
save(all_abun.table, file = "06.Species_genetic_association/RData/NSCLC_RCC_all_abun.table.RData")

all_abun.table$Species<-France_info$Short_name[match(all_abun.table$Species, France_info$organism)]

sv_assoc_id<-paste(all_adonis.table$Species, all_adonis.table$Prog,sep = "_")
abun_assoc_id<-paste(all_abun.table$Species, all_abun.table$Prog, sep = "_")

all_abun.table<-all_abun.table[match(sv_assoc_id,abun_assoc_id),]

colnames(all_abun.table)<-paste("Abun",colnames(all_abun.table),sep = ".")
save(all_abun.table, file = "06.Species_genetic_association/RData/NSCLC_RCC_all_abun.table.RData")

species.table<-cbind(all_adonis.table,all_abun.table)

species.table$sv.MetaSigAssoc<-rep('No', nrow(species.table))
species.table$sv.MetaSigAssoc[species.table$p< 0.05]<-'Yes'
#species.table$sv.MetaSigAssoc[species.table$fdr< 0.05]<-'Yes'
species.table$sv.MetaSigAssoc[is.na(species.table$p)==T]<-'Unknown'

species.table$abun.MetaSigAssoc<-rep('No', nrow(species.table))
species.table$abun.MetaSigAssoc[species.table$Abun.p< 0.05]<-'Yes'
#species.table$abun.MetaSigAssoc[species.table$Abun.fdr< 0.05]<-'Yes'

species.sig.table<-species.table[,] %>%
  .[.$sv.MetaSigAssoc=="Yes" | .$abun.MetaSigAssoc == "Yes",]

write.table(species.sig.table, "06.Species_genetic_association/NSCLC_RCC_species.sig.tsv",sep = "\t", 
   col.names = T, row.names = F, quote = F)

write.csv(species.sig.table, "06.Species_genetic_association/NSCLC_RCC_species.sig.csv",row.names = F)

unique(species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]$Species)
unique(species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]$Prog)

unique(species.sig.table[species.sig.table$abun.MetaSigAssoc=="Yes",]$Species)
unique(species.sig.table[species.sig.table$abun.MetaSigAssoc=="Yes",]$Prog)

species.sv.sig.table<-species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]


### 3.4 Venn diagram
species_count<-table(species.table$sv.MetaSigAssoc,species.table$abun.MetaSigAssoc)

#pdf("06.Species_genetic_association/species_Prog_count.venn.pdf", width = 2, height = 2)
tiff(file = "pics/FS2C_NSCLC_RCC_species_count.venn.tiff", width =600, height =600, res =300) 

draw.pairwise.venn(species_count[3,1]+species_count[3,2],
                   species_count[1,2]+species_count[2,2]+species_count[3,2],
                   species_count[3,2], 
                   category = c("Genetics", "Abundance"), lty = rep("blank",2), 
                   fill =c("#4472c4", "#00b050"), alpha = rep(0.5, 2), cat.pos = c(0, 0), 
                   cat.cex = c(0.5, 0.5),cat.dist = rep(0.025, 2), scaled = F)

dev.off()
 

##############################################################################
###  3.5 heatmap
load("06.Species_genetic_association/RData/NSCLC_RCC_all_adonis.table.Rdata")
load("06.Species_genetic_association/RData/NSCLC_RCC_all_abun.table.RData")

## matrix
adonis_res.P<-all_adonis.table[,c("Prog", "Species", "p")] %>%
  spread("Prog", "p")
adonis_res.P<-data.frame(adonis_res.P, row.names = "Species")

all_abun_res.P<-all_abun.table[,c("Abun.Prog", "Abun.Species", "Abun.p")] %>%
  spread("Abun.Prog", "Abun.p")

all_abun_res.P<-data.frame(all_abun_res.P, row.names = "Abun.Species")
all_abun_res.P<-all_abun_res.P[match(rownames(adonis_res.P),rownames(all_abun_res.P)),
                                                     match(colnames(adonis_res.P),colnames(all_abun_res.P))]

# color
species.color<-matrix(0, nrow = nrow(adonis_res.P), ncol = ncol(adonis_res.P))
rownames(species.color)<-rownames(adonis_res.P)
colnames(species.color)<-colnames(adonis_res.P)

species.color[adonis_res.P<0.05&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 1        #"purple", both assoc
species.color[adonis_res.P>=0.05&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 2       #"yellow", abundance assoc
species.color[adonis_res.P<0.05&all_abun_res.P>=0.05&!is.na(all_abun_res.P)]<- 3       #"blue", genetic assoc
species.color[is.na(adonis_res.P)==T&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 4   #"grey", none assoc


### association number 
species.asso_number<-matrix(0, nrow = nrow(adonis_res.P), ncol = ncol(adonis_res.P))
rownames(species.asso_number)<-rownames(adonis_res.P)
colnames(species.asso_number)<-colnames(adonis_res.P)

species.asso_number[adonis_res.P<0.05&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 2        #both assoc
species.asso_number[adonis_res.P>=0.05&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 1       #abundance assoc
species.asso_number[adonis_res.P<0.05&all_abun_res.P>=0.05&!is.na(all_abun_res.P)]<- 1       #genetic assoc
species.asso_number[is.na(adonis_res.P)==T&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 0   #none assoc

species.asso<-apply(species.asso_number,1,sum)
species.asso<-species.asso[which(species.asso>0)]
species.asso<-species.asso[order(species.asso,na.last = TRUE, decreasing = TRUE)]
species_order<-names(species.asso)

## get plot tables
##species.sig.table<-species.table[species.table$fdr<0.05 | species.table$Abun.fdr<0.05,]
species.sig.table<-species.table[species.table$p<0.05 | species.table$Abun.p<0.05,]
species.sig.table<-subset(species.sig.table,is.na(sv.MetaSigAssoc)==F,drop=T)

##write.csv(species.color,"test.csv")
species.plot<- species.sig.table$Prog %>%
  as.character(.) %>%
  na.omit %>%
  .[!duplicated(.)]

species.plot.spe <- species.sig.table$Species %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

species.color.plot    <- species.color %>%
  .[match(species.plot.spe,rownames(.)),match(species.plot,colnames(.))]

colnames(species.color.plot)<-c("Response (NSCLC)","Response (RCC)")
species.color.plot<-species.color.plot[species_order,]

tiff(file = "pics/F4A_NSCLC_RCC_specie_prognosis.heatmap_N50.tiff", width =1000, height =2800, res =300) 

heatmap.2(species.color.plot, 
          col=colorRampPalette(c("#e9e9e9","#993399","#E18727FF", "#0072B5FF"))(4), #
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",
          cexCol =1.5, srtCol = 62, cexRow =1.5,
          labRow=as.expression(lapply(rownames(species.color.plot), function(a) bquote(italic(.(a))))),
          colsep = c(1:(ncol(species.color.plot)-1)), rowsep = c(1:(nrow(species.color.plot)-1)),
          sepcolor="#30303a",sepwidth = c(0.02,0.02),
          key = F,
          lmat=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.1, 3.5),lwid=c(3, 7, 4.2 ),key.title = NA,
          margins=c(16,7))

dev.off()

#mycolor6
#hist(mtcars$mpg,
#     breaks = 15,
#     col = "#666666",
#)


#######################################################################################
### NSCLC

species_short_name<-"R.lactaris"
species_dis<-PRJNA751792_msv_dist_std[[paste("msv_",France_info$organism[match(species_short_name,France_info$Short_name)],sep = "")]]

prog_name<-"response_code"
clin_PRJNA751792<-subset(all_basic,dataset=="PRJNA751792",drop=T)
prog_vec<-clin_PRJNA751792[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Non-responder (SD/PD)"
prog_vec_input[prog_vec_input==1]<-"Responder (CR/PR)"

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ggtitle("NSCLC (DerosaL_2022):\nR.lactaris\nadonis R2=0.02;P-value:0.025")+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#BC3C29FF","#6F99ADFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 10,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))


#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F4C_R.lactaris_NSCLC_response.tiff", width =1800, height =2300, res =600) 
print(p_spe_pcoa)
dev.off()


##########################################################################################
## some species 
species_short_name<-"B.uniformis"
species_dis<-PRJEB22863_msv_dist_std[[paste("msv_",France_info$organism[match(species_short_name,France_info$Short_name)],sep = "")]]

prog_name<-"response_code"
clin_PRJEB22863<-subset(all_basic,dataset=="PRJEB22863",drop=T)
clin_PRJEB22863_NSCLC<-subset(all_basic,cancer_type=="NSCLC",drop=T)
prog_vec<-clin_PRJEB22863_NSCLC[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Non-responser"
prog_vec_input[prog_vec_input==1]<-"Responser"

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ggtitle("NSCLC (DerosaL_2022):\nR.lactaris\nadonis R2=0.02;P-value:0.003")+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#BC3C29FF","#6F99ADFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 10,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))


#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F4C_R.lactaris_NSCLC_response.tiff", width =1800, height =2300, res =600) 
print(p_spe_pcoa)
dev.off()



##########################################################################################
## some species 
species_short_name<-"B.adolescentis"
species_dis<-PRJNA751792_msv_dist_std[[paste("msv_",France_info$organism[match(species_short_name,France_info$Short_name)],sep = "")]]

prog_name<-"response_code"
clin_PRJNA751792<-subset(all_basic,dataset=="PRJNA751792",drop=T)

prog_vec<-clin_PRJNA751792[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Non-responser"
prog_vec_input[prog_vec_input==1]<-"Responser"

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ggtitle("NSCLC (DerosaL_2022):\nB.adolescentis\nadonis R2=0.05;P-value:0.014")+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#BC3C29FF","#6F99ADFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 10,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))


#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F4C_B.adolescentis_NSCLC_response.tiff", width =1800, height =2300, res =600) 
print(p_spe_pcoa)
dev.off()






##########################################################################################
## some species 
species_short_name<-"A.muciniphila"
species_dis<-PRJEB22863_msv_dist_std[[paste("msv_",France_info$organism[match(species_short_name,France_info$Short_name)],sep = "")]]

prog_name<-"response_code"
clin_PRJEB22863<-subset(all_basic,dataset=="PRJEB22863",drop=T)

clin_PRJEB22863_RCC<-subset(clin_PRJEB22863,cancer_type=="RCC",drop=T)
prog_vec<-clin_PRJEB22863_RCC[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Non-responder (SD/PD)"
prog_vec_input[prog_vec_input==1]<-"Responder (CR/PR)"

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ggtitle("RCC (RoutyB_2018):\nA.muciniphila\nadonis R2=0.07;P-value:0.046")+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#BC3C29FF","#6F99ADFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 10,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))


#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F4C_A.muciniphila_RCC_response.tiff", width =1800, height =2300, res =600) 
print(p_spe_pcoa)
dev.off()



############################################################################################
##  输出表格 (Mel)

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")
USA_UK_info<- read.table("01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R/06.Species_genetic_association")

species.sig<-read.csv("Mel_species.sig.csv",header=T)

Melanoma_irAEs<-subset(species.sig,Prog=="Melanoma_irAEs",drop=T)
Melanoma_os_medians<-subset(species.sig,Prog=="Melanoma_os_medians",drop=T)
Melanoma_pfs_12_months<-subset(species.sig,Prog=="Melanoma_pfs_12_months",drop=T)
Melanoma_response<-subset(species.sig,Prog=="Melanoma_response",drop=T)

################# Melanoma_irAEs
mel_irAEs_adonis<-read.csv("mel_irAEs_adonis.csv",header=T)
mel_irAEs_abun<-read.csv("mel_irAEs_abun.csv",header=T)
colnames(mel_irAEs_adonis)[1]<-"Species"
mel_irAEs_abun$Species<-USA_UK_info$Short_name[match(mel_irAEs_abun$McCullochJA_2022.X, USA_UK_info$organism)]

Melanoma_irAEs_all<-merge(Melanoma_irAEs,mel_irAEs_adonis,by.x="Species",by.y="Species")
Melanoma_irAEs_all<-merge(Melanoma_irAEs_all,mel_irAEs_abun,by.x="Species",by.y="Species")

write.csv(Melanoma_irAEs_all,"table/Melanoma_irAEs_all.csv",row.names=F)


##################### Melanoma_pfs_12_months
mel_pfs_12_months_adonis<-read.csv("mel_pfs_12_months_adonis.csv",header=T)
mel_pfs_12_months_abun<-read.csv("mel_pfs_12_months_abun.csv",header=T)
colnames(mel_pfs_12_months_adonis)[1]<-"Species"
mel_pfs_12_months_abun$Species<-USA_UK_info$Short_name[match(mel_pfs_12_months_abun$McCullochJA_2022.X, USA_UK_info$organism)]

Melanoma_pfs_12_months_all<-merge(Melanoma_pfs_12_months,mel_pfs_12_months_adonis,by.x="Species",by.y="Species")
Melanoma_pfs_12_months_all<-merge(Melanoma_pfs_12_months_all,mel_pfs_12_months_abun,by.x="Species",by.y="Species")

write.csv(Melanoma_pfs_12_months_all,"table/Melanoma_pfs_12_months_all.csv",row.names=F)


##################### Melanoma_response
mel_response_adonis<-read.csv("mel_resp_adonis.csv",header=T)
mel_response_abun<-read.csv("mel_resp_abun.csv",header=T)
colnames(mel_response_adonis)[1]<-"Species"
mel_response_abun$Species<-USA_UK_info$Short_name[match(mel_response_abun$McCullochJA_2022.X,USA_UK_info$organism)]

Melanoma_response_all<-merge(Melanoma_response,mel_response_adonis,by.x="Species",by.y="Species")
Melanoma_response_all<-merge(Melanoma_response_all,mel_response_abun,by.x="Species",by.y="Species")

write.csv(Melanoma_response_all,"table/Melanoma_response_all.csv",row.names=F)



############################################################################################
##  输出表格 (NSCLC_RCC)

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")
France_info<- read.table("01.cleanData/SV_info/France_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R/06.Species_genetic_association")

species.sig<-read.csv("NSCLC_RCC_species.sig.csv",header=T)

nsclc_os<-subset(species.sig,Prog=="NSCLC_os_medians",drop=T)
nsclc_response<-subset(species.sig,Prog=="NSCLC_response",drop=T)
rcc_response<-subset(species.sig,Prog=="RCC_response",drop=T)


######################### NSCLC response 
nsclc_response_adonis<-read.csv("nsclc_resp_adonis.csv",header=T)
nsclc_response_abun<-read.csv("nsclc_resp_abun.csv",header=T)
colnames(nsclc_response_adonis)[1]<-"Species"
nsclc_response_abun$Species<-France_info$Short_name[match(nsclc_response_abun$DerosaL_2022.X, France_info$organism)]

nsclc_response_all<-merge(nsclc_response,nsclc_response_adonis,by.x="Species",by.y="Species")
nsclc_response_all<-merge(nsclc_response_all,nsclc_response_abun,by.x="Species",by.y="Species")

write.csv(nsclc_response_all,"table/nsclc_response_all.csv",row.names=F)


######################### rcc response 
rcc_response_adonis<-read.csv("rcc_resp_adonis.csv",header=T)
rcc_response_abun<-read.csv("rcc_resp_abun.csv",header=T)
colnames(rcc_response_adonis)[1]<-"Species"
rcc_response_abun$Species<-France_info$Short_name[match(rcc_response_abun$RoutyB_2018_RCC.X, France_info$organism)]

rcc_response_all<-merge(rcc_response,rcc_response_adonis,by.x="Species",by.y="Species")
rcc_response_all<-merge(rcc_response_all,rcc_response_abun,by.x="Species",by.y="Species")

write.csv(rcc_response_all,"table/rcc_response_all.csv",row.names=F)




#######################################################
##  calculate number

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R/06.Species_genetic_association")

Mel_species.sig<-read.csv("Mel_species.sig.csv",header=T)
NSCLC_RCC_species.sig<-read.csv("NSCLC_RCC_species.sig.csv",header=T)

mel_abun<-subset(Mel_species.sig,abun.MetaSigAssoc=="Yes",drop=T)
nsclc_abun<-subset(NSCLC_RCC_species.sig,abun.MetaSigAssoc=="Yes",drop=T)

abun<-rbind(mel_abun,nsclc_abun)

nrow(abun)
length(unique(abun$Species))

##write.csv(abun,"abun_test.csv")



mel_sv<-subset(Mel_species.sig,sv.MetaSigAssoc=="Yes",drop=T)
nsclc_sv<-subset(NSCLC_RCC_species.sig,sv.MetaSigAssoc=="Yes",drop=T)

sv<-rbind(mel_sv,nsclc_sv)

nrow(sv)
length(unique(sv$Species))

##write.csv(sv,"SV_test.csv")



