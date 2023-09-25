### Microbial GWA 
### 2023-9-20
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
options(scipen = 200)

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)
library("R.utils")
library(showtext)
showtext_auto()

#library("survival")
#library("ggplot2")
#library("survminer")
#library(gridExtra)
#library("grid")
library(reshape2)	  
#library("RColorBrewer")
#library("plyr")



########################################################
#### 3.2 Replication 
#### 3.2.1 vSV

### Correlation of effect size between different cancer types 
###  response 

RCC_dsv_resp<-read.csv("07.Microbial_GWAS/datasets/rcc_dsv_resp_pvalue_adjust.csv",header=T)
RCC_vsv_resp<-read.csv("07.Microbial_GWAS/datasets/rcc_vsv_resp_pvalue_adjust.csv",header=T)
RCC_resp<-rbind(RCC_dsv_resp,RCC_vsv_resp)
RCC_resp_sub<-RCC_resp[,c("X","OR")]
colnames(RCC_resp_sub)<-c("id","RCC_OR")

NSCLC_dsv_resp<-read.csv("07.Microbial_GWAS/datasets/nsclc_dsv_resp_meta_pvalue.csv",header=T)
NSCLC_vsv_resp<-read.csv("07.Microbial_GWAS/datasets/nsclc_vsv_resp_meta_pvalue.csv",header=T)
NSCLC_resp<-rbind(NSCLC_dsv_resp,NSCLC_vsv_resp)
NSCLC_resp_sub<-NSCLC_resp[,c("X","effect_meta")]
NSCLC_resp_sub$NSCLC_OR<-exp(NSCLC_resp$effect_meta)
colnames(NSCLC_resp_sub)<-c("id","NSCLC_OR")

melanoma_dsv_resp<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_resp_meta_pvalue.csv",header=T)
melanoma_vsv_resp<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_resp_meta_pvalue.csv",header=T)
melanoma_resp<-rbind(melanoma_dsv_resp,melanoma_vsv_resp)
melanoma_resp_sub<-melanoma_resp[,c("X","effect_meta")]
melanoma_resp_sub$melanoma_OR<-exp(melanoma_resp$effect_meta)
colnames(melanoma_resp_sub)<-c("id","melanoma_OR")

resp<-merge(RCC_resp_sub,NSCLC_resp_sub,by.x="id",by.y="id")
resp<-merge(resp,melanoma_resp_sub,by.x="id",by.y="id")

####################################################################
## NSCLC and melanoma 

resp<-subset(resp,melanoma_OR<30,drop=T)
resp<-subset(resp,NSCLC_OR<30,drop=T)
es_cor<-cor.test(resp$NSCLC_OR, resp$melanoma_OR,method = "spearman")

text_r<-paste("R=",round(es_cor$estimate,digits = 2),
              ", P=0.80",
              sep = "")

p_resp_mel<-ggplot(resp) + 
  geom_point(aes(NSCLC_OR, melanoma_OR),alpha = 0.4) +
  annotate("text",x = c(2), y = c(7.5), label = c(text_r),hjust = 0,size=8)+
  geom_hline(yintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  geom_vline(xintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  xlab('OR (NSCLC)')+
  ylab('OR (melanoma)')+
  scale_color_viridis()+
 theme_bw(base_size = 17)+
  theme(legend.position = 'none',axis.text=element_text(size=20))


####################################################################
## Melanoma and RCC

resp<-subset(resp,RCC_OR<30,drop=T)
resp<-subset(resp,melanoma_OR<30,drop=T)
es_cor<-cor.test(resp$melanoma_OR, resp$RCC_OR,method = "spearman")

text_r<-paste("R=",round(es_cor$estimate,digits = 2),
              ", P=0.02",
              sep = "")

p_resp_rcc<-ggplot(resp) + 
  geom_point(aes(melanoma_OR, RCC_OR),alpha = 0.4) +
  annotate("text",x = c(2), y = c(15), label = c(text_r),hjust = 0,size=8)+
  geom_hline(yintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  geom_vline(xintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  xlab('OR (Melanoma)')+
  ylab('OR (RCC)')+
  scale_color_viridis()+
 theme_bw(base_size = 17)+
  theme(legend.position = 'none',axis.text=element_text(size=20))


####################################################################
## RCC and NSCLC

resp<-subset(resp,NSCLC_OR<30,drop=T)
resp<-subset(resp,RCC_OR<30,drop=T)
es_cor<-cor.test(resp$RCC_OR, resp$NSCLC_OR,method = "spearman")

text_r<-paste("R=",round(es_cor$estimate,digits = 2),
              ", P=0.82",
              sep = "")

p_resp_NSCLC<-ggplot(resp) + 
  geom_point(aes(RCC_OR, NSCLC_OR),alpha = 0.4) +
  annotate("text",x = c(5), y = c(5), label = c(text_r),hjust = 0,size=8)+
  geom_hline(yintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  geom_vline(xintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  xlab('OR (RCC)')+
  ylab('OR (NSCLC)')+
  scale_color_viridis()+
 theme_bw(base_size = 17)+
  theme(legend.position = 'none',axis.text=element_text(size=20))


#######################################################################################
## os 

NSCLC_dsv_os<-read.csv("07.Microbial_GWAS/datasets/nsclc_dsv_os_pvalue_adjust.csv",header=T)
NSCLC_vsv_os<-read.csv("07.Microbial_GWAS/datasets/nsclc_vsv_os_pvalue_adjust.csv",header=T)
NSCLC_os<-rbind(NSCLC_dsv_os,NSCLC_vsv_os)
NSCLC_os_sub<-NSCLC_os[,c("X","HR")]
colnames(NSCLC_os_sub)<-c("id","NSCLC_HR")

melanoma_dsv_os<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_os_pvalue_adjust.csv",header=T)
melanoma_vsv_os<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_os_pvalue_adjust.csv",header=T)
melanoma_os<-rbind(melanoma_dsv_os,melanoma_vsv_os)
melanoma_os_sub<-melanoma_os[,c("X","HR")]
colnames(melanoma_os_sub)<-c("id","melanoma_HR")

os<-merge(NSCLC_os_sub,melanoma_os_sub,by.x="id",by.y="id")

####################################################################
## overall and melanoma 

os<-subset(os,melanoma_HR<30,drop=T)
os<-subset(os,NSCLC_HR<30,drop=T)
es_cor<-cor.test(os$NSCLC_HR, os$melanoma_HR,method = "spearman")

text_r<-paste("R=",round(es_cor$estimate,digits = 2),
              ", P=0.37",
              sep = "")

p_os_melanoma<-ggplot(os) + 
  geom_point(aes(NSCLC_HR, melanoma_HR),alpha = 0.4) +
  annotate("text",x = c(1), y = c(11), label = c(text_r),hjust = 0,size=8)+
  geom_hline(yintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  geom_vline(xintercept=0, linetype="dashed",  color = "red", size=1, alpha = 0.5)+
  xlab('HR (NSCLC)')+
  ylab('HR (melanoma)')+
  scale_color_viridis()+
  theme_bw(base_size = 17)+
  theme(legend.position = 'none',axis.text=element_text(size=20))



#### 3.2.3 Combine figure
## plot

p_title_replic <- ggdraw() + 
    #draw_label(
    #   'Replication',
    #  fontface = 'bold', x = 0, hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))

p_sv_prog_replic<-plot_grid(
    p_title_replic, 
    plot_grid(p_resp_mel,p_resp_rcc,p_resp_NSCLC,p_os_melanoma,
              rel_widths = c(1, 1, 1,1),align = 'hv',
              labels = c("SVs-response OR", "SVs-response OR","SVs-response OR",
                         "SVs-os HR "),
              ncol = 2,label_size	=18,vjust = 0),
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 2)
  )



tiff(file = "pics/FS6_replicate.tiff", width =1500, height =1500, res =300) 
print(p_sv_prog_replic)
dev.off()




###################### 1.2 Inputs
USA_UK_info<- read.table("01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

USA_UK_vsv_info<-read.table("01.cleanData/SV_info/USA_UK_vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
USA_UK_dsv_info<-read.table("01.cleanData/SV_info/USA_UK_dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")


###################################################################################################
### 3 Results
### 3.1 Clean results
#### 3.1.1 vSV associations

## vsv
## merge result tables


###################### 1.2 Inputs
France_info<- read.table("01.cleanData/SV_info/France_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)
USA_UK_info<- read.table("01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

France_vsv_info<-read.table("01.cleanData/SV_info/France_vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
France_dsv_info<-read.table("01.cleanData/SV_info/France_dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

USA_UK_vsv_info<-read.table("01.cleanData/SV_info/USA_UK_vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
USA_UK_dsv_info<-read.table("01.cleanData/SV_info/USA_UK_dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")



mel_vsv_resp_meta_pvalue<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_resp_meta_pvalue_adjust.csv",header=T)
mel_vsv_resp_meta_pvalue_sig<-subset(mel_vsv_resp_meta_pvalue,pvalue_meta<0.05,drop=T)
melanoma_vsv_resp_final<-left_join(mel_vsv_resp_meta_pvalue_sig, USA_UK_vsv_info, by = c("X" = "SV_Name"))

melanoma_vsv_resp_select<-read.csv("07.Microbial_GWAS/melanoma_vsv_response_adjAbun.sig.anno.csv",header=T)
melanoma_vsv_resp_select$Figure<-"T"
melanoma_vsv_resp_select<-melanoma_vsv_resp_select[,c("X","Figure")]
melanoma_vsv_resp_final<-merge(melanoma_vsv_resp_final,melanoma_vsv_resp_select,by.x="X",by.y="X",all.x=T)

write.csv(melanoma_vsv_resp_final,"07.Microbial_GWAS/combine/melanoma_vsv_resp.csv",row.names=F)


mel_vsv_pfs_12_months_meta_pvalue<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_pfs_12_months_meta_pvalue_adjust.csv",header=T)
mel_vsv_pfs_12_months_meta_pvalue_sig<-subset(mel_vsv_pfs_12_months_meta_pvalue,pvalue_meta<0.05,drop=T)
melanoma_vsv_pfs_12_months_final<-left_join(mel_vsv_pfs_12_months_meta_pvalue_sig, USA_UK_vsv_info, by = c("X" = "SV_Name"))

melanoma_vsv_pfs_select<-read.csv("07.Microbial_GWAS/melanoma_vsv_pfs_12_months_adjAbun.sig.anno.csv",header=T)
melanoma_vsv_pfs_select$Figure<-"T"
melanoma_vsv_pfs_select<-melanoma_vsv_pfs_select[,c("X","Figure")]
melanoma_vsv_pfs_12_months_final<-merge(melanoma_vsv_pfs_12_months_final,melanoma_vsv_pfs_select,by.x="X",by.y="X",all.x=T)

write.csv(melanoma_vsv_pfs_12_months_final,"07.Microbial_GWAS/combine/melanoma_vsv_pfs_12_months.csv",row.names=F)


mel_vsv_irAEs_pvalue_adjust<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_irAEs_pvalue_adjust.csv",header=T)
mel_vsv_irAEs_pvalue_adjust_sig<-subset(mel_vsv_irAEs_pvalue_adjust,p<0.05,drop=T)
melanoma_vsv_irAEs_final<-left_join(mel_vsv_irAEs_pvalue_adjust_sig, USA_UK_vsv_info, by = c("X" = "SV_Name"))
melanoma_vsv_irAEs_final$Figure<-"NA"
write.csv(melanoma_vsv_irAEs_final,"07.Microbial_GWAS/combine/melanoma_vsv_irAEs.csv",row.names=F)


mel_vsv_os_pvalue_adjust<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_os_pvalue_adjust.csv",header=T)
mel_vsv_os_pvalue_adjust_sig<-subset(mel_vsv_os_pvalue_adjust,p<0.05,drop=T)
melanoma_vsv_os_final<-left_join(mel_vsv_os_pvalue_adjust_sig, USA_UK_vsv_info, by = c("X" = "SV_Name"))
melanoma_vsv_os_final$Figure<-"NA"
write.csv(melanoma_vsv_os_final,"07.Microbial_GWAS/combine/melanoma_vsv_os.csv",row.names=F)


mel_dsv_resp_meta_pvalue<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_resp_meta_pvalue_adjust.csv",header=T)
mel_dsv_resp_meta_pvalue_sig<-subset(mel_dsv_resp_meta_pvalue,pvalue_meta<0.05,drop=T)
melanoma_dsv_resp_final<-left_join(mel_dsv_resp_meta_pvalue_sig, USA_UK_dsv_info, by = c("X" = "SV_Name"))

melanoma_dsv_resp_select<-read.csv("07.Microbial_GWAS/melanoma_dsv_response_adjAbun.sig.anno.csv",header=T)
melanoma_dsv_resp_select$Figure<-"T"
melanoma_dsv_resp_select<-melanoma_dsv_resp_select[,c("X","Figure")]
melanoma_dsv_resp_final<-merge(melanoma_dsv_resp_final,melanoma_dsv_resp_select,by.x="X",by.y="X",all.x=T)

write.csv(melanoma_dsv_resp_final,"07.Microbial_GWAS/combine/melanoma_dsv_resp.csv",row.names=F)


mel_dsv_pfs_12_months_meta_pvalue<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_pfs_12_months_meta_pvalue_adjust.csv",header=T)
mel_dsv_pfs_12_months_meta_pvalue_sig<-subset(mel_dsv_pfs_12_months_meta_pvalue,pvalue_meta<0.05,drop=T)
melanoma_dsv_pfs_12_months_final<-left_join(mel_dsv_pfs_12_months_meta_pvalue_sig, USA_UK_dsv_info, by = c("X" = "SV_Name"))

melanoma_dsv_pfs_select<-read.csv("07.Microbial_GWAS/melanoma_dsv_pfs_12_months_adjAbun.sig.anno.csv",header=T)
melanoma_dsv_pfs_select$Figure<-"T"
melanoma_dsv_pfs_select<-melanoma_dsv_pfs_select[,c("X","Figure")]
melanoma_dsv_pfs_12_months_final<-merge(melanoma_dsv_pfs_12_months_final,melanoma_dsv_pfs_select,by.x="X",by.y="X",all.x=T)

write.csv(melanoma_dsv_pfs_12_months_final,"07.Microbial_GWAS/combine/melanoma_dsv_pfs_12_months.csv",row.names=F)


mel_dsv_irAEs_pvalue_adjust<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_irAEs_pvalue_adjust.csv",header=T)
mel_dsv_irAEs_pvalue_adjust_sig<-subset(mel_dsv_irAEs_pvalue_adjust,p<0.05,drop=T)
melanoma_dsv_irAEs_final<-left_join(mel_dsv_irAEs_pvalue_adjust_sig, USA_UK_dsv_info, by = c("X" = "SV_Name"))
melanoma_dsv_irAEs_final$Figure<-"NA"

write.csv(melanoma_dsv_irAEs_final,"07.Microbial_GWAS/combine/melanoma_dsv_irAEs.csv",row.names=F)


mel_dsv_os_pvalue_adjust<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_os_pvalue_adjust.csv",header=T)
mel_dsv_os_pvalue_adjust_sig<-subset(mel_dsv_os_pvalue_adjust,p<0.05,drop=T)
melanoma_dsv_os_final<-left_join(mel_dsv_os_pvalue_adjust_sig, USA_UK_dsv_info, by = c("X" = "SV_Name"))
melanoma_dsv_os_final$Figure<-"NA"

write.csv(melanoma_dsv_os_final,"07.Microbial_GWAS/combine/melanoma_dsv_os.csv",row.names=F)


#### the species count

dsv_os<-read.csv("07.Microbial_GWAS/combine/melanoma_dsv_os.csv",header=T)$Taxonomy_Name
vsv_os<-read.csv("07.Microbial_GWAS/combine/melanoma_vsv_os.csv",header=T)$Taxonomy_Name
os<-c(dsv_os,vsv_os)
length(unique(os))

vsv_resp<-read.csv("07.Microbial_GWAS/combine/melanoma_vsv_resp.csv",header=T)$Taxonomy_Name
dsv_resp<-read.csv("07.Microbial_GWAS/combine/melanoma_dsv_resp.csv",header=T)$Taxonomy_Name
resp<-c(vsv_resp,dsv_resp)
length(unique(resp))



#######################  NSCLC

###################### 1.2 Inputs
France_info<- read.table("01.cleanData/SV_info/France_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

France_vsv_info<-read.table("01.cleanData/SV_info/France_vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
France_dsv_info<-read.table("01.cleanData/SV_info/France_dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

nsclc_vsv_resp_meta_pvalue<-read.csv("07.Microbial_GWAS/datasets/nsclc_vsv_resp_meta_pvalue_adjust.csv",header=T)
nsclc_vsv_resp_meta_pvalue_sig<-subset(nsclc_vsv_resp_meta_pvalue,pvalue_meta<0.05,drop=T)
nsclc_vsv_resp_final<-left_join(nsclc_vsv_resp_meta_pvalue_sig, France_vsv_info, by = c("X" = "SV_Name"))

nsclc_vsv_resp_select<-read.csv("07.Microbial_GWAS/NSCLC_vsv_response_adjAbun.sig.anno.csv",header=T)
nsclc_vsv_resp_select$Figure<-"T"
nsclc_vsv_resp_select<-nsclc_vsv_resp_select[,c("X","Figure")]
nsclc_vsv_resp_final<-merge(nsclc_vsv_resp_final,nsclc_vsv_resp_select,by.x="X",by.y="X",all.x=T)

write.csv(nsclc_vsv_resp_final,"07.Microbial_GWAS/combine/nsclc_vsv_resp.csv",row.names=F)

nsclc_vsv_os_pvalue_adjust<-read.csv("07.Microbial_GWAS/datasets/nsclc_vsv_os_pvalue_adjust.csv",header=T)
nsclc_vsv_os_pvalue_adjust_sig<-subset(nsclc_vsv_os_pvalue_adjust,p<0.05,drop=T)
nsclc_vsv_os_final<-left_join(nsclc_vsv_os_pvalue_adjust_sig, France_vsv_info, by = c("X" = "SV_Name"))

nsclc_vsv_os_select<-read.csv("07.Microbial_GWAS/NSCLC_vsv_os_adjAbun.sig.anno.csv",header=T)
nsclc_vsv_os_select$Figure<-"T"
nsclc_vsv_os_select<-nsclc_vsv_os_select[,c("X","Figure")]
nsclc_vsv_os_final<-merge(nsclc_vsv_os_final,nsclc_vsv_os_select,by.x="X",by.y="X",all.x=T)

write.csv(nsclc_vsv_os_final,"07.Microbial_GWAS/combine/nsclc_vsv_os.csv",row.names=F)


nsclc_dsv_resp_meta_pvalue<-read.csv("07.Microbial_GWAS/datasets/nsclc_dsv_resp_meta_pvalue_adjust.csv",header=T)
nsclc_dsv_resp_meta_pvalue_sig<-subset(nsclc_dsv_resp_meta_pvalue,pvalue_meta<0.05,drop=T)
nsclc_dsv_resp_final<-left_join(nsclc_dsv_resp_meta_pvalue_sig, France_dsv_info, by = c("X" = "SV_Name"))

nsclc_dsv_resp_select<-read.csv("07.Microbial_GWAS/NSCLC_dsv_response_adjAbun.sig.anno.csv",header=T)
nsclc_dsv_resp_select$Figure<-"T"
nsclc_dsv_resp_select<-nsclc_dsv_resp_select[,c("X","Figure")]
nsclc_dsv_resp_final<-merge(nsclc_dsv_resp_final,nsclc_dsv_resp_select,by.x="X",by.y="X",all.x=T)

write.csv(nsclc_dsv_resp_final,"07.Microbial_GWAS/combine/nsclc_dsv_resp.csv",row.names=F)

nsclc_dsv_os_pvalue_adjust<-read.csv("07.Microbial_GWAS/datasets/nsclc_dsv_os_pvalue_adjust.csv",header=T)
nsclc_dsv_os_pvalue_adjust_sig<-subset(nsclc_dsv_os_pvalue_adjust,p<0.05,drop=T)
nsclc_dsv_os_final<-left_join(nsclc_dsv_os_pvalue_adjust_sig, France_dsv_info, by = c("X" = "SV_Name"))

nsclc_dsv_os_select<-read.csv("07.Microbial_GWAS/NSCLC_dsv_os_adjAbun.sig.anno.csv",header=T)
nsclc_dsv_os_select$Figure<-"T"
nsclc_dsv_os_select<-nsclc_dsv_os_select[,c("X","Figure")]
nsclc_dsv_os_final<-merge(nsclc_dsv_os_final,nsclc_dsv_os_select,by.x="X",by.y="X",all.x=T)

write.csv(nsclc_dsv_os_final,"07.Microbial_GWAS/combine/nsclc_dsv_os.csv",row.names=F)



#######################  RCC
rcc_vsv_resp_meta_pvalue<-read.csv("07.Microbial_GWAS/datasets/rcc_vsv_resp_pvalue_adjust.csv",header=T)
rcc_vsv_resp_meta_pvalue_sig<-subset(rcc_vsv_resp_meta_pvalue,p<0.05,drop=T)
rcc_vsv_resp_final<-left_join(rcc_vsv_resp_meta_pvalue_sig, France_vsv_info, by = c("X" = "SV_Name"))

rcc_vsv_resp_final$Figure<-"NA"
write.csv(rcc_vsv_resp_final,"07.Microbial_GWAS/combine/rcc_vsv_resp.csv",row.names=F)


#######################  RCC
rcc_dsv_resp_meta_pvalue<-read.csv("07.Microbial_GWAS/datasets/rcc_dsv_resp_pvalue_adjust.csv",header=T)
rcc_dsv_resp_meta_pvalue_sig<-subset(rcc_dsv_resp_meta_pvalue,p<0.05,drop=T)
rcc_dsv_resp_final<-left_join(rcc_dsv_resp_meta_pvalue_sig, France_dsv_info, by = c("X" = "SV_Name"))

rcc_dsv_resp_select<-read.csv("07.Microbial_GWAS/RCC_dsv_response_adjAbun.sig.anno.csv",header=T)
rcc_dsv_resp_select$Figure<-"T"
rcc_dsv_resp_select<-rcc_dsv_resp_select[,c("X","Figure")]
rcc_dsv_resp_final<-merge(rcc_dsv_resp_final,rcc_dsv_resp_select,by.x="X",by.y="X",all.x=T)

write.csv(rcc_dsv_resp_final,"07.Microbial_GWAS/combine/rcc_dsv_resp.csv",row.names=F)



