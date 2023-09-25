### Microbial GWA 
### 2023-9-19
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
#options(scipen = 200)

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)
#library("R.utils")
#library(showtext)
#showtext_auto()

#library("survival")
#library("ggplot2")
#library("survminer")
#library(gridExtra)
#library("grid")
library(reshape2)	  
#library("RColorBrewer")
#library("plyr")
library("metafor")


###################### 1.2 Inputs
France_info<- read.table("01.cleanData/SV_info/France_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)
USA_UK_info<- read.table("01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

France_vsv_info<-read.table("01.cleanData/SV_info/France_vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
France_dsv_info<-read.table("01.cleanData/SV_info/France_dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

USA_UK_vsv_info<-read.table("01.cleanData/SV_info/USA_UK_vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
USA_UK_dsv_info<-read.table("01.cleanData/SV_info/USA_UK_dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

France_basic<-read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")
USA_UK_basic<-read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv")

all_basic<-rbind(France_basic,USA_UK_basic)

PRJEB22863_basic<-subset(all_basic,dataset=="PRJEB22863",drop=T)
PRJNA397906_basic<-subset(all_basic,dataset=="PRJNA397906",drop=T)
PRJNA751792_basic<-subset(all_basic,dataset=="PRJNA751792",drop=T)
PRJNA762360_basic<-subset(all_basic,dataset=="PRJNA762360",drop=T)
PRJNA770295_basic<-subset(all_basic,dataset=="PRJNA770295",drop=T)
PRJEB43119_basic<-subset(all_basic,dataset=="PRJEB43119",drop=T)
PRJNA541981_basic<-subset(all_basic,dataset=="PRJNA541981",drop=T)
PRJEB22863_NSCLC_basic<-subset(PRJEB22863_basic,cancer_type=="NSCLC",drop=T)
PRJEB22863_RCC_basic<-subset(PRJEB22863_basic,cancer_type=="RCC",drop=T)


load("01.cleanData/SV_all/dsgv_PRJEB22863.RData")
load("01.cleanData/SV_all/vsgv_PRJEB22863.RData")

PRJEB22863_NSCLC_dsgv<-PRJEB22863_dsgv[row.names(PRJEB22863_NSCLC_basic),]
PRJEB22863_NSCLC_vsgv<-PRJEB22863_vsgv[row.names(PRJEB22863_NSCLC_basic),]

PRJEB22863_RCC_dsgv<-PRJEB22863_dsgv[row.names(PRJEB22863_RCC_basic),]
PRJEB22863_RCC_vsgv<-PRJEB22863_vsgv[row.names(PRJEB22863_RCC_basic),]

load("01.cleanData/SV_all/dsgv_PRJNA762360.RData")
load("01.cleanData/SV_all/vsgv_PRJNA762360.RData")

load("01.cleanData/SV_all/dsgv_PRJNA397906.RData")
load("01.cleanData/SV_all/vsgv_PRJNA397906.RData")

load("01.cleanData/SV_all/dsgv_PRJNA751792.RData")
load("01.cleanData/SV_all/vsgv_PRJNA751792.RData")

load("01.cleanData/SV_all/dsgv_PRJNA770295.RData")
load("01.cleanData/SV_all/vsgv_PRJNA770295.RData")

load("01.cleanData/SV_all/dsgv_PRJEB43119.RData")
load("01.cleanData/SV_all/vsgv_PRJEB43119.RData")

load("01.cleanData/SV_all/dsgv_PRJNA541981.RData")
load("01.cleanData/SV_all/vsgv_PRJNA541981.RData")


PRJEB22863_abun<-read.table("01.cleanData/mbio_all/PRJEB22863_SV_species_abun.tsv", check.names = F) 
PRJEB22863_NSCLC_abun<-PRJEB22863_abun[row.names(PRJEB22863_NSCLC_basic),]
PRJEB22863_RCC_abun<-PRJEB22863_abun[row.names(PRJEB22863_RCC_basic),]
PRJNA397906_abun<-read.table("01.cleanData/mbio_all/PRJNA397906_SV_species_abun.tsv", check.names = F)  
PRJNA751792_abun<-read.table("01.cleanData/mbio_all/PRJNA751792_SV_species_abun.tsv", check.names = F) 
PRJNA762360_abun<-read.table("01.cleanData/mbio_all/PRJNA762360_SV_species_abun.tsv", check.names = F) 
PRJNA770295_abun<-read.table("01.cleanData/mbio_all/PRJNA770295_SV_species_abun.tsv", check.names = F)
PRJEB43119_abun<-read.table("01.cleanData/mbio_all/PRJEB43119_SV_species_abun.tsv", check.names = F)
PRJNA541981_abun<-read.table("01.cleanData/mbio_all/PRJNA541981_SV_species_abun.tsv", check.names = F)


## 2 Associations between SVs and prognosis
### 2.1 Preparation
if (!dir.exists("07.Microbial_GWAS")) {dir.create("07.Microbial_GWAS")}
if (!dir.exists("07.Microbial_GWAS/datasets")) {dir.create("07.Microbial_GWAS/datasets")}


# PRJEB22863
PRJEB22863_abun$Others<-1-rowSums(PRJEB22863_abun)
PRJEB22863_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB22863_abun)), transform="clr") %>%as.data.frame
PRJEB22863_abun_clr <- PRJEB22863_abun_clr[match(rownames(PRJEB22863_abun), rownames(PRJEB22863_abun_clr)),]
rownames(PRJEB22863_abun_clr) <- rownames(PRJEB22863_abun)
PRJEB22863_covar<-cbind(PRJEB22863_basic,PRJEB22863_abun_clr)

# PRJEB22863 NSCLC
PRJEB22863_NSCLC_abun$Others<-1-rowSums(PRJEB22863_NSCLC_abun)
PRJEB22863_NSCLC_abun_clr<-abundances(x=as.data.frame(na.omit(PRJEB22863_NSCLC_abun)), transform="clr") %>%as.data.frame
PRJEB22863_NSCLC_abun_clr <- PRJEB22863_NSCLC_abun_clr[match(rownames(PRJEB22863_NSCLC_abun), rownames(PRJEB22863_NSCLC_abun_clr)),]
rownames(PRJEB22863_NSCLC_abun_clr) <- rownames(PRJEB22863_NSCLC_abun)
PRJEB22863_NSCLC_covar<-cbind(PRJEB22863_NSCLC_basic,PRJEB22863_NSCLC_abun_clr)

# PRJEB22863 RCC
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


############################################################################################################################
###  for dataset (PRJNA762360)
###  Logistic regression model 1
#Linear model with PRJNA762360 covariates: age, gender.

PRJNA762360_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600_2"<-PRJNA762360_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600"

PRJNA762360_covar <- c("age_bins","gender_num")
PRJNA762360_prog<-PRJNA762360_basic[,c("response_code","irAEs","pfs_12_months")]

PRJNA762360_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJNA762360_prog,PRJNA762360_vsgv,PRJNA762360_basic,PRJNA762360_covar,USA_UK_info,"response_code")
PRJNA762360_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJNA762360_prog,PRJNA762360_dsgv,PRJNA762360_basic,PRJNA762360_covar,USA_UK_info,"response_code")

write.csv(PRJNA762360_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA762360_vsv_resp_res.csv",row.names=F)
write.csv(PRJNA762360_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA762360_dsv_resp_res.csv",row.names=F)

PRJNA762360_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJNA762360_prog,PRJNA762360_vsgv,PRJNA762360_basic,PRJNA762360_covar,USA_UK_info,"irAEs")
PRJNA762360_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJNA762360_prog,PRJNA762360_dsgv,PRJNA762360_basic,PRJNA762360_covar,USA_UK_info,"irAEs")

write.csv(PRJNA762360_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA762360_vsv_irAEs_res.csv",row.names=F)
write.csv(PRJNA762360_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA762360_dsv_irAEs_res.csv",row.names=F)

PRJNA762360_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJNA762360_prog,PRJNA762360_vsgv,PRJNA762360_basic,PRJNA762360_covar,USA_UK_info,"pfs_12_months")
PRJNA762360_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJNA762360_prog,PRJNA762360_dsgv,PRJNA762360_basic,PRJNA762360_covar,USA_UK_info,"pfs_12_months")

write.csv(PRJNA762360_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA762360_vsv_pfs_12_months_res.csv",row.names=F)
write.csv(PRJNA762360_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA762360_dsv_pfs_12_months_res.csv",row.names=F)


### Cox regression model 1
#Survival model with PRJNA762360_covariates: age, gender

PRJNA762360_surv<-PRJNA762360_basic[,c("os","os_event")]
PRJNA762360_vsv_os_cox_res<-cox_btw_mats_noadjAbun_vsv(PRJNA762360_surv,PRJNA762360_vsgv,PRJNA762360_basic,PRJNA762360_covar,USA_UK_info,"os","os_event")
PRJNA762360_dsv_os_cox_res<-cox_btw_mats_noadjAbun_dsv(PRJNA762360_surv,PRJNA762360_dsgv,PRJNA762360_basic,PRJNA762360_covar,USA_UK_info,"os","os_event")

write.csv(PRJNA762360_dsv_os_cox_res, file = "07.Microbial_GWAS/datasets/PRJNA762360_dsv_os_res.csv",row.names=F)
write.csv(PRJNA762360_vsv_os_cox_res, file = "07.Microbial_GWAS/datasets/PRJNA762360_vsv_os_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJNA397906)
###  Logistic regression model 1
#Linear model with PRJNA397906 covariates: age, gender
# 1 indicate the organism name column of the info file

PRJNA397906_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600_2"<-PRJNA397906_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600"

PRJNA397906_covar <- c("age_bins","gender_num")
PRJNA397906_prog<-PRJNA397906_basic[,c("response_code","irAEs")]

PRJNA397906_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJNA397906_prog,PRJNA397906_vsgv,PRJNA397906_basic,PRJNA397906_covar,USA_UK_info,"response_code")
PRJNA397906_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJNA397906_prog,PRJNA397906_dsgv,PRJNA397906_basic,PRJNA397906_covar,USA_UK_info,"response_code")

write.csv(PRJNA397906_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA397906_vsv_resp_res.csv",row.names=F)
write.csv(PRJNA397906_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA397906_dsv_resp_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJNA770295)
###  Logistic regression model 1
#Linear model with PRJNA770295 covariates: age, gender.
# 1 indicate the organism name column of the info file

PRJNA770295_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600_2"<-PRJNA770295_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600"

PRJNA770295_covar <- c("age_bins","gender_num")
PRJNA770295_prog<-PRJNA770295_basic[,c("response_code","pfs_12_months")]

PRJNA770295_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJNA770295_prog,PRJNA770295_vsgv,PRJNA770295_basic,PRJNA770295_covar,USA_UK_info,"response_code")
PRJNA770295_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJNA770295_prog,PRJNA770295_dsgv,PRJNA770295_basic,PRJNA770295_covar,USA_UK_info,"response_code")

write.csv(PRJNA770295_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA770295_vsv_resp_res.csv",row.names=F)
write.csv(PRJNA770295_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA770295_dsv_resp_res.csv",row.names=F)

PRJNA770295_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJNA770295_prog,PRJNA770295_vsgv,PRJNA770295_basic,PRJNA770295_covar,USA_UK_info,"pfs_12_months")
PRJNA770295_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJNA770295_prog,PRJNA770295_dsgv,PRJNA770295_basic,PRJNA770295_covar,USA_UK_info,"pfs_12_months")

write.csv(PRJNA770295_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA770295_vsv_pfs_12_months_res.csv",row.names=F)
write.csv(PRJNA770295_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA770295_dsv_pfs_12_months_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJEB43119)
###  Logistic regression model 1
#Linear model with PRJEB43119 covariates: age, gender
# 1 indicate the organism name column of the info file

PRJEB43119_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600_2"<-PRJEB43119_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600"

PRJEB43119_covar <- c("age_bins","gender_num")
PRJEB43119_prog<-PRJEB43119_basic[,c("response_code","pfs_12_months")]

PRJEB43119_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJEB43119_prog,PRJEB43119_vsgv,PRJEB43119_basic,PRJEB43119_covar,USA_UK_info,"response_code")
PRJEB43119_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJEB43119_prog,PRJEB43119_dsgv,PRJEB43119_basic,PRJEB43119_covar,USA_UK_info,"response_code")

write.csv(PRJEB43119_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB43119_vsv_resp_res.csv",row.names=F)
write.csv(PRJEB43119_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB43119_dsv_resp_res.csv",row.names=F)

PRJEB43119_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJEB43119_prog,PRJEB43119_vsgv,PRJEB43119_basic,PRJEB43119_covar,USA_UK_info,"pfs_12_months")
PRJEB43119_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJEB43119_prog,PRJEB43119_dsgv,PRJEB43119_basic,PRJEB43119_covar,USA_UK_info,"pfs_12_months")

write.csv(PRJEB43119_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB43119_vsv_pfs_12_months_res.csv",row.names=F)
write.csv(PRJEB43119_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB43119_dsv_pfs_12_months_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJNA541981)
###  Logistic regression model 1
#Linear model with PRJNA541981 covariates: age, gender.
# 1 indicate the organism name column of the info file

PRJNA541981_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600_2"<-PRJNA541981_vsgv$"Phascolarctobacterium sp. CAG:207:1597_1600"

PRJNA541981_covar <- c("age_bins","gender_num")
PRJNA541981_prog<-PRJNA541981_basic[,c("pfs_12_months","pfs_12_months")]

PRJNA541981_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJNA541981_prog,PRJNA541981_vsgv,PRJNA541981_basic,PRJNA541981_covar,USA_UK_info,"pfs_12_months")
PRJNA541981_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJNA541981_prog,PRJNA541981_dsgv,PRJNA541981_basic,PRJNA541981_covar,USA_UK_info,"pfs_12_months")

write.csv(PRJNA541981_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA541981_vsv_pfs_12_months_res.csv",row.names=F)
write.csv(PRJNA541981_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA541981_dsv_pfs_12_months_res.csv",row.names=F)


### Cox regression model 1
#Survival model with PRJNA541981_covariates: age, gender

PRJNA541981_surv<-PRJNA541981_basic[,c("os","os_event")]
PRJNA541981_vsv_os_cox_res<-cox_btw_mats_noadjAbun_vsv(PRJNA541981_surv,PRJNA541981_vsgv,PRJNA541981_basic,PRJNA541981_covar,USA_UK_info,"os","os_event")
PRJNA541981_dsv_os_cox_res<-cox_btw_mats_noadjAbun_dsv(PRJNA541981_surv,PRJNA541981_dsgv,PRJNA541981_basic,PRJNA541981_covar,USA_UK_info,"os","os_event")

write.csv(PRJNA541981_dsv_os_cox_res, file = "07.Microbial_GWAS/datasets/PRJNA541981_dsv_os_res.csv",row.names=F)
write.csv(PRJNA541981_vsv_os_cox_res, file = "07.Microbial_GWAS/datasets/PRJNA541981_vsv_os_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJNA751792)
###  Logistic regression model 1
#Linear model with PRJNA751792 covariates: age, gender.
# 1 indicate the organism name column of the info file

PRJNA751792_covar <- c("age_bins","gender_num")
PRJNA751792_prog<-PRJNA751792_basic[,c("response_code","response_code")]

PRJNA751792_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJNA751792_prog,PRJNA751792_vsgv,PRJNA751792_basic,PRJNA751792_covar,France_info,"response_code")
PRJNA751792_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJNA751792_prog,PRJNA751792_dsgv,PRJNA751792_basic,PRJNA751792_covar,France_info,"response_code")

write.csv(PRJNA751792_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA751792_vsv_resp_res.csv",row.names=F)
write.csv(PRJNA751792_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJNA751792_dsv_resp_res.csv",row.names=F)

### Cox regression model 1
#Survival model with PRJNA751792_covariates: age, gender, BMI, read count and corresponding species relative abundance.

PRJNA751792_surv<-PRJNA751792_basic[,c("os","os_event")]
PRJNA751792_vsv_os_cox_res<-cox_btw_mats_noadjAbun_vsv(PRJNA751792_surv,PRJNA751792_vsgv,PRJNA751792_basic,PRJNA751792_covar,France_info,"os","os_event")
PRJNA751792_dsv_os_cox_res<-cox_btw_mats_noadjAbun_vsv(PRJNA751792_surv,PRJNA751792_dsgv,PRJNA751792_basic,PRJNA751792_covar,France_info,"os","os_event")

write.csv(PRJNA751792_dsv_os_cox_res, file = "07.Microbial_GWAS/datasets/PRJNA751792_dsv_os_res.csv",row.names=F)
write.csv(PRJNA751792_vsv_os_cox_res, file = "07.Microbial_GWAS/datasets/PRJNA751792_vsv_os_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJEB22863)
###  Logistic regression model 1
#Linear model with PRJEB22863 covariates: age, gender.
# 1 indicate the organism name column of the info file

PRJEB22863_covar <- c("age_bins","gender_num")
PRJEB22863_prog<-PRJEB22863_basic[,c("response_code","irAEs")]

PRJEB22863_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJEB22863_prog,PRJEB22863_vsgv,PRJEB22863_basic,PRJEB22863_covar,France_info,"response_code")
PRJEB22863_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJEB22863_prog,PRJEB22863_dsgv,PRJEB22863_basic,PRJEB22863_covar,France_info,"response_code")

write.csv(PRJEB22863_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB22863_vsv_resp_res.csv",row.names=F)
write.csv(PRJEB22863_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB22863_dsv_resp_res.csv",row.names=F)

############################################################################################################################
###  for dataset (PRJEB22863) NSCLC
###  Logistic regression model 1
#Linear model with PRJEB22863 covariates: age, gender.
# 1 indicate the organism name column of the info file

PRJEB22863_covar <- c("age_bins","gender_num")
PRJEB22863_prog<-PRJEB22863_NSCLC_basic[,c("response_code","irAEs")]

PRJEB22863_NSCLC_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJEB22863_prog,PRJEB22863_NSCLC_vsgv,PRJEB22863_NSCLC_basic,PRJEB22863_covar,France_info,"response_code")
PRJEB22863_NSCLC_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJEB22863_prog,PRJEB22863_NSCLC_dsgv,PRJEB22863_NSCLC_basic,PRJEB22863_covar,France_info,"response_code")

write.csv(PRJEB22863_NSCLC_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB22863_NSCLC_vsv_resp_res.csv",row.names=F)
write.csv(PRJEB22863_NSCLC_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB22863_NSCLC_dsv_resp_res.csv",row.names=F)


############################################################################################################################
###  for dataset (PRJEB22863) RCC
###  Logistic regression model 1
#Linear model with PRJEB22863 covariates: age, gender.
# 1 indicate the organism name column of the info file

PRJEB22863_covar <- c("age_bins","gender_num")
PRJEB22863_prog<-PRJEB22863_RCC_basic[,c("response_code","irAEs")]

PRJEB22863_RCC_vsv_lg_res<-logistic_btw_mats_noadjAbun_vsv(PRJEB22863_prog,PRJEB22863_RCC_vsgv,PRJEB22863_RCC_basic,PRJEB22863_covar,France_info,"response_code")
PRJEB22863_RCC_dsv_lg_res<-logistic_btw_mats_noadjAbun_dsv(PRJEB22863_prog,PRJEB22863_RCC_dsgv,PRJEB22863_RCC_basic,PRJEB22863_covar,France_info,"response_code")

write.csv(PRJEB22863_RCC_vsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB22863_RCC_vsv_resp_res.csv",row.names=F)
write.csv(PRJEB22863_RCC_dsv_lg_res, file = "07.Microbial_GWAS/datasets/PRJEB22863_RCC_dsv_resp_res.csv",row.names=F)


#################################################################################################################################
### dsv meta-analysis

################ response 

PRJEB22863_dsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB22863_dsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJEB22863_dsv_resp_res.table)[-1]<-c(paste("RoutyB_2018.",colnames(PRJEB22863_dsv_resp_res.table)[2:5],sep = ""))

PRJEB22863_NSCLC_dsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB22863_NSCLC_dsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJEB22863_NSCLC_dsv_resp_res.table)[-1]<-c(paste("RoutyB_2018.",colnames(PRJEB22863_NSCLC_dsv_resp_res.table)[2:5],sep = ""))

PRJNA397906_dsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA397906_dsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA397906_dsv_resp_res.table)[-1]<-c(paste("FrankelAE_2017.",colnames(PRJNA397906_dsv_resp_res.table)[2:5],sep = ""))

PRJNA751792_dsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA751792_dsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA751792_dsv_resp_res.table)[-1]<-c(paste("DerosaL_2022.",colnames(PRJNA751792_dsv_resp_res.table)[2:5],sep = ""))

PRJNA762360_dsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA762360_dsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA762360_dsv_resp_res.table)[-1]<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_dsv_resp_res.table)[2:5],sep = ""))

PRJNA770295_dsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA770295_dsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA770295_dsv_resp_res.table)[-1]<-c(paste("SpencerCN_2021.",colnames(PRJNA770295_dsv_resp_res.table)[2:5],sep = ""))

PRJEB43119_dsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB43119_dsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJEB43119_dsv_resp_res.table)[-1]<-c(paste("LeeKA_2022.",colnames(PRJEB43119_dsv_resp_res.table)[2:5],sep = ""))

### mel
mel_dsv_resp_res.table<-merge(PRJNA397906_dsv_resp_res.table,PRJNA762360_dsv_resp_res.table,by.x="X",by.y="X")
mel_dsv_resp_res.table<-merge(mel_dsv_resp_res.table,PRJNA770295_dsv_resp_res.table,by.x="X",by.y="X")
mel_dsv_resp_res.table<-merge(mel_dsv_resp_res.table,PRJEB43119_dsv_resp_res.table,by.x="X",by.y="X")

combine_meta_p_lr(mel_dsv_resp_res.table,"07.Microbial_GWAS/datasets/mel_dsv_resp_meta_pvalue.csv")
mel_dsv_resp<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_resp_meta_pvalue.csv",header=T)
multi_adjust(mel_dsv_resp,"07.Microbial_GWAS/datasets/mel_dsv_resp_meta_pvalue_adjust.csv")

### nsclc
nsclc_dsv_resp_res.table<-merge(PRJNA751792_dsv_resp_res.table,PRJEB22863_NSCLC_dsv_resp_res.table,by.x="X",by.y="X")
combine_meta_p_lr(nsclc_dsv_resp_res.table,"07.Microbial_GWAS/datasets/nsclc_dsv_resp_meta_pvalue.csv")

nsclc_dsv_resp<-read.csv("07.Microbial_GWAS/datasets/nsclc_dsv_resp_meta_pvalue.csv",header=T)
multi_adjust(nsclc_dsv_resp,"07.Microbial_GWAS/datasets/nsclc_dsv_resp_meta_pvalue_adjust.csv")

### rcc (only one datasets)
PRJEB22863_RCC_dsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB22863_RCC_dsv_resp_res.csv",header=T)
write.csv(PRJEB22863_RCC_dsv_resp_res.table,"07.Microbial_GWAS/datasets/rcc_dsv_resp_pvalue_adjust.csv")


#####################################################
## pfs 12 months 

PRJNA762360_dsv_pfs_12_months_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA762360_dsv_pfs_12_months_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA762360_dsv_pfs_12_months_res.table)[-1]<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_dsv_pfs_12_months_res.table)[2:5],sep = ""))

PRJNA770295_dsv_pfs_12_months_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA770295_dsv_pfs_12_months_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA770295_dsv_pfs_12_months_res.table)[-1]<-c(paste("SpencerCN_2021.",colnames(PRJNA770295_dsv_pfs_12_months_res.table)[2:5],sep = ""))

PRJEB43119_dsv_pfs_12_months_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB43119_dsv_pfs_12_months_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJEB43119_dsv_pfs_12_months_res.table)[-1]<-c(paste("LeeKA_2022.",colnames(PRJEB43119_dsv_pfs_12_months_res.table)[2:5],sep = ""))

PRJNA541981_dsv_pfs_12_months_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA541981_dsv_pfs_12_months_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA541981_dsv_pfs_12_months_res.table)[-1]<-c(paste("PetersBA_2019.",colnames(PRJNA541981_dsv_pfs_12_months_res.table)[2:5],sep = ""))

### mel
mel_dsv_pfs_12_months_res.table<-merge(PRJNA762360_dsv_pfs_12_months_res.table,PRJEB43119_dsv_pfs_12_months_res.table,by.x="X",by.y="X")
mel_dsv_pfs_12_months_res.table<-merge(mel_dsv_pfs_12_months_res.table,PRJNA770295_dsv_pfs_12_months_res.table,by.x="X",by.y="X")
mel_dsv_pfs_12_months_res.table<-merge(mel_dsv_pfs_12_months_res.table,PRJNA541981_dsv_pfs_12_months_res.table,by.x="X",by.y="X")

combine_meta_p_lr(mel_dsv_pfs_12_months_res.table,"07.Microbial_GWAS/datasets/mel_dsv_pfs_12_months_meta_pvalue.csv")
mel_dsv_pfs_12_months<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_pfs_12_months_meta_pvalue.csv",header=T)
multi_adjust(mel_dsv_pfs_12_months,"07.Microbial_GWAS/datasets/mel_dsv_pfs_12_months_meta_pvalue_adjust.csv")


###############################################################################################
### os 
### mel (only one dataset)
PRJNA762360_dsv_os_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA762360_dsv_os_res.csv",header=T)
write.csv(PRJNA762360_dsv_os_res.table,"07.Microbial_GWAS/datasets/mel_dsv_os_pvalue_adjust.csv",row.names=F)

### nsclc (only one dataset)
PRJNA751792_dsv_os_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA751792_dsv_os_res.csv",header=T)
write.csv(PRJNA751792_dsv_os_res.table,"07.Microbial_GWAS/datasets/nsclc_dsv_os_pvalue_adjust.csv",row.names=F)

###################################
## irAEs
#### mel (only one dataset)
PRJNA762360_dsv_irAEs_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA762360_dsv_irAEs_res.csv",header=T)
write.csv(PRJNA762360_dsv_irAEs_res.table,"07.Microbial_GWAS/datasets/mel_dsv_irAEs_pvalue_adjust.csv")


###############################################################################
### vsv meta-analysis

PRJEB22863_vsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB22863_vsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJEB22863_vsv_resp_res.table)[-1]<-c(paste("RoutyB_2018.",colnames(PRJEB22863_vsv_resp_res.table)[2:5],sep = ""))

PRJEB22863_NSCLC_vsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB22863_NSCLC_vsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJEB22863_NSCLC_vsv_resp_res.table)[-1]<-c(paste("RoutyB_2018.",colnames(PRJEB22863_NSCLC_vsv_resp_res.table)[2:5],sep = ""))

PRJNA397906_vsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA397906_vsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA397906_vsv_resp_res.table)[-1]<-c(paste("FrankelAE_2017.",colnames(PRJNA397906_vsv_resp_res.table)[2:5],sep = ""))

PRJNA751792_vsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA751792_vsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA751792_vsv_resp_res.table)[-1]<-c(paste("DerosaL_2022.",colnames(PRJNA751792_vsv_resp_res.table)[2:5],sep = ""))

PRJNA762360_vsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA762360_vsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA762360_vsv_resp_res.table)[-1]<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_vsv_resp_res.table)[2:5],sep = ""))

PRJNA770295_vsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA770295_vsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA770295_vsv_resp_res.table)[-1]<-c(paste("SpencerCN_2021.",colnames(PRJNA770295_vsv_resp_res.table)[2:5],sep = ""))

PRJEB43119_vsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB43119_vsv_resp_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJEB43119_vsv_resp_res.table)[-1]<-c(paste("LeeKA_2022.",colnames(PRJEB43119_vsv_resp_res.table)[2:5],sep = ""))


###   p meta  (all samples)
#combine_meta_p_lr(cbind_vsv_resp_res.table,"07.Microbial_GWAS/datasets/all_vsv_resp_meta_pvalue.csv")
#multi_adjust("07.Microbial_GWAS/datasets/all_vsv_resp_meta_pvalue.csv","07.Microbial_GWAS/datasets/all_vsv_resp_meta_pvalue_adjust.csv")
#write.csv(cbind_vsv_resp_res.table,"07.Microbial_GWAS/datasets/all_vsv_resp_pvalue.csv",row.names=F)

### mel
mel_vsv_resp_res.table<-merge(PRJNA397906_vsv_resp_res.table,PRJNA762360_vsv_resp_res.table,by.x="X",by.y="X")
mel_vsv_resp_res.table<-merge(mel_vsv_resp_res.table,PRJNA770295_vsv_resp_res.table,by.x="X",by.y="X")
mel_vsv_resp_res.table<-merge(mel_vsv_resp_res.table,PRJEB43119_vsv_resp_res.table,by.x="X",by.y="X")

combine_meta_p_lr(mel_vsv_resp_res.table,"07.Microbial_GWAS/datasets/mel_vsv_resp_meta_pvalue.csv")
mel_vsv_resp<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_resp_meta_pvalue.csv",header=T)
multi_adjust(mel_vsv_resp,"07.Microbial_GWAS/datasets/mel_vsv_resp_meta_pvalue_adjust.csv")

### nsclc
nsclc_vsv_resp_res.table<-merge(PRJNA751792_vsv_resp_res.table,PRJEB22863_NSCLC_vsv_resp_res.table,by.x="X",by.y="X")
#write.csv(nsclc_vsv_resp_res.table,"07.Microbial_GWAS/datasets/nsclc_vsv_resp_pvalue.csv",row.names=F)
combine_meta_p_lr(nsclc_vsv_resp_res.table,"07.Microbial_GWAS/datasets/nsclc_vsv_resp_meta_pvalue.csv")

nsclc_vsv_resp<-read.csv("07.Microbial_GWAS/datasets/nsclc_vsv_resp_meta_pvalue.csv",header=T)
multi_adjust(nsclc_vsv_resp,"07.Microbial_GWAS/datasets/nsclc_vsv_resp_meta_pvalue_adjust.csv")

### rcc (only one datasets)
PRJEB22863_RCC_vsv_resp_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB22863_RCC_vsv_resp_res.csv",header=T)
write.csv(PRJEB22863_RCC_vsv_resp_res.table,"07.Microbial_GWAS/datasets/rcc_vsv_resp_pvalue_adjust.csv")


#####################################################
## pfs 12 months 

PRJNA762360_vsv_pfs_12_months_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA762360_vsv_pfs_12_months_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA762360_vsv_pfs_12_months_res.table)[-1]<-c(paste("McCullochJA_2022.",colnames(PRJNA762360_vsv_pfs_12_months_res.table)[2:5],sep = ""))

PRJNA770295_vsv_pfs_12_months_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA770295_vsv_pfs_12_months_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA770295_vsv_pfs_12_months_res.table)[-1]<-c(paste("SpencerCN_2021.",colnames(PRJNA770295_vsv_pfs_12_months_res.table)[2:5],sep = ""))

PRJEB43119_vsv_pfs_12_months_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJEB43119_vsv_pfs_12_months_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJEB43119_vsv_pfs_12_months_res.table)[-1]<-c(paste("LeeKA_2022.",colnames(PRJEB43119_vsv_pfs_12_months_res.table)[2:5],sep = ""))

PRJNA541981_vsv_pfs_12_months_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA541981_vsv_pfs_12_months_res.csv",header=T)[,c(seq(1,4),7)]
colnames(PRJNA541981_vsv_pfs_12_months_res.table)[-1]<-c(paste("PetersBA_2019.",colnames(PRJNA541981_vsv_pfs_12_months_res.table)[2:5],sep = ""))

### mel
mel_vsv_pfs_12_months_res.table<-merge(PRJNA762360_vsv_pfs_12_months_res.table,PRJEB43119_vsv_pfs_12_months_res.table,by.x="X",by.y="X")
mel_vsv_pfs_12_months_res.table<-merge(mel_vsv_pfs_12_months_res.table,PRJNA770295_vsv_pfs_12_months_res.table,by.x="X",by.y="X")
mel_vsv_pfs_12_months_res.table<-merge(mel_vsv_pfs_12_months_res.table,PRJNA541981_vsv_pfs_12_months_res.table,by.x="X",by.y="X")

combine_meta_p_lr(mel_vsv_pfs_12_months_res.table,"07.Microbial_GWAS/datasets/mel_vsv_pfs_12_months_meta_pvalue.csv")
mel_vsv_pfs_12_months<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_pfs_12_months_meta_pvalue.csv",header=T)
multi_adjust(mel_vsv_pfs_12_months,"07.Microbial_GWAS/datasets/mel_vsv_pfs_12_months_meta_pvalue_adjust.csv")

###############################################################################################
### os 
### mel 
PRJNA762360_vsv_os_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA762360_vsv_os_res.csv",header=T)
write.csv(PRJNA762360_vsv_os_res.table,"07.Microbial_GWAS/datasets/mel_vsv_os_pvalue_adjust.csv",row.names=F)

### nsclc (only one dataset)
PRJNA751792_vsv_os_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA751792_vsv_os_res.csv",header=T)
write.csv(PRJNA751792_vsv_os_res.table,"07.Microbial_GWAS/datasets/nsclc_vsv_os_pvalue_adjust.csv",row.names=F)

###################################
## irAEs
#### mel (only one dataset)
PRJNA762360_vsv_irAEs_res.table<-read.csv("07.Microbial_GWAS/datasets/PRJNA762360_vsv_irAEs_res.csv",header=T)
write.csv(PRJNA762360_vsv_irAEs_res.table,"07.Microbial_GWAS/datasets/mel_vsv_irAEs_pvalue_adjust.csv",row.names=F)




