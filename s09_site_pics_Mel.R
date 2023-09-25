### Survival and barplot for sites
### 2023-9-18
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")

## 1 Preparation
##rm(list = ls())

library("survival")
library("survminer")
library(gridExtra)
library("grid")
#library(reshape2)	  
#library("RColorBrewer")
library("plyr")

#install.packages("patchwork")
library("patchwork")
source("functions.R")

###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/USA_UK_vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/USA_UK_dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv",header=T)
mel_basic<-subset(all_basic,cancer_type=="melanoma",drop=T)

load("01.cleanData/SV_all/USA_UK_dsgv.RData")
load("01.cleanData/SV_all/USA_UK_vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub

#######################################################################################
### clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Eggerthella lenta DSM 2243:2572_2575;2575_2578",
"Faecalibacterium cf. prausnitzii KLE1255:458_459",
"Subdoligranulum sp. 4_3_54A2FAA:2121_2122",
"Akkermansia muciniphila ATCC BAA-835:863_865"
)]


### draw pictures
### response
#Akkermansia muciniphila ATCC BAA-835:863_865
#Subdoligranulum sp. 4_3_54A2FAA:2121_2122

### pfs_12_months
## Akkermansia muciniphila ATCC BAA-835:863_865
##Faecalibacterium cf. prausnitzii KLE1255:458_459

########################################################################################
##  response (melano)
dsgv_pic$id<-row.names(dsgv_pic)
pic_dsv<-merge(mel_basic,dsgv_pic,by.x="id",by.y="id")

#########################################################################
#####  "Subdoligranulum sp. 4_3_54A2FAA:2121_2122"

pic_dsv$temp<-pic_dsv$"Subdoligranulum sp. 4_3_54A2FAA:2121_2122"
lg_res_model <- glm(response_code~temp+age,family=binomial(logit), data=pic_dsv)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(response_code)==F,drop=T)
table(pic_sub$study)

FrankelAE_2017_pic<-subset(pic_sub,study=="FrankelAE_2017",drop=T)
LeeKA_2022_pic<-subset(pic_sub,study=="LeeKA_2022",drop=T)
McCullochJA_2022_pic<-subset(pic_sub,study=="McCullochJA_2022",drop=T)
SpencerCN_2021_pic<-subset(pic_sub,study=="SpencerCN_2021",drop=T)

FrankelAE_2017_ratio<-dsv_resp(FrankelAE_2017_pic,"FrankelAE_2017")
LeeKA_2022_ratio<-dsv_resp(LeeKA_2022_pic,"LeeKA_2022")
McCullochJA_2022_ratio<-dsv_resp(McCullochJA_2022_pic,"McCullochJA_2022")
SpencerCN_2021_ratio<-dsv_resp(SpencerCN_2021_pic,"SpencerCN_2021")

###par(family = "serif")
table(FrankelAE_2017_pic$temp)

p1<-ggplot(FrankelAE_2017_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("FrankelAE_\n2017")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=10)","Non-Delection (n=14)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=11,face="plain")
   )

table(McCullochJA_2022_pic$temp)
p2<-ggplot(McCullochJA_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("McCullochJA_\n2022")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=15)","Non-Delection (n=11)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

table(SpencerCN_2021_pic$temp)
p3<-ggplot(SpencerCN_2021_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("SpencerCN_\n2021")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=16)","Non-Delection (n=17)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

table(LeeKA_2022_pic$temp)
p4<-ggplot(LeeKA_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("LeeKA_\n2022")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y ='')+
   geom_text(aes(y=new_col, label=value),size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=16)","Non-Delection (n=27)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),axis.text.y = element_blank(),
   plot.title=element_text(size=11,face="plain")
   )

p_whole<-p1+p2+p3+p4+plot_layout(ncol=4) + plot_annotation(title = "Melenoma\ndSV Subdoligranulum sp. 4_3_54A2FAA:2121_2122 kbp\nMeta OR=0.28;Meta p=2.56e-3")

#p_whole<-plot_grid(p1,p2,p3,p4,
#              rel_widths = c(0.55,0.45,0.45,0.7),align = 'hv',
#              ncol = 4,label_size= 10,vjust = 0)

tiff(file = "pics/sites/F5_dsv_mel_resp_Subdoligranulumsp4_3_54A2FAA_2121_2122_barplot.tiff", width =1900, height =1300, res =300) 
p_whole
dev.off()


########################################################################################
##  response Akkermansia muciniphila ATCC BAA-835:863_865

pic_dsv$temp<-pic_dsv$"Akkermansia muciniphila ATCC BAA-835:863_865"
lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic_dsv)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic_dsv,is.na(temp)==F,drop=T)
table(pic_sub$study)

FrankelAE_2017_pic<-subset(pic_sub,study=="FrankelAE_2017",drop=T)
LeeKA_2022_pic<-subset(pic_sub,study=="LeeKA_2022",drop=T)
McCullochJA_2022_pic<-subset(pic_sub,study=="McCullochJA_2022",drop=T)
SpencerCN_2021_pic<-subset(pic_sub,study=="SpencerCN_2021",drop=T)

FrankelAE_2017_ratio<-dsv_resp(FrankelAE_2017_pic,"FrankelAE_2017")
LeeKA_2022_ratio<-dsv_resp(LeeKA_2022_pic,"LeeKA_2022")
McCullochJA_2022_ratio<-dsv_resp(McCullochJA_2022_pic,"McCullochJA_2022")
SpencerCN_2021_ratio<-dsv_resp(SpencerCN_2021_pic,"SpencerCN_2021")

table(McCullochJA_2022_pic$temp)
p1<-ggplot(McCullochJA_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("McCullochJA_\n2022")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=3)","Non-Delection (n=14)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=11,face="plain")
   )

table(SpencerCN_2021_pic$temp)
p2<-ggplot(SpencerCN_2021_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("SpencerCN_\n2021")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=20)","Non-Delection (n=10)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size =10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

table(LeeKA_2022_pic$temp)
p3<-ggplot(LeeKA_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("LeeKA_\n2022")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y ='')+
   geom_text(aes(y=new_col, label=value),size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=27)","Non-Delection (n=29)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),axis.text.y = element_blank(),
   plot.title=element_text(size=11,face="plain")
   )

p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "Melenoma\ndSV Akkermansia muciniphila ATCC BAA-835:863~865 kbp\nMeta OR=0.24;Meta p=2.73e-3")

#p_whole<-plot_grid(p1,p2,p3,p4,
#              rel_widths = c(0.55,0.45,0.45,0.7),align = 'hv',
#              ncol = 4,label_size= 10,vjust = 0)

tiff(file = "pics/sites/F5_dsv_mel_resp_A.muciniphila863_865_barplot.tiff", width =1700, height =1300, res =300) 
p_whole
dev.off()


###################################################################################
### pfs_12_months
## Akkermansia muciniphila ATCC BAA-835:863_865

pic_dsv$temp<-pic_dsv$"Akkermansia muciniphila ATCC BAA-835:863_865"
lg_res_model <- glm(pfs_12_months~temp,family=binomial(logit), data=pic_dsv)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(pfs_12_months)==F,drop=T)
table(pic_sub$study)

LeeKA_2022_pic<-subset(pic_sub,study=="LeeKA_2022",drop=T)
McCullochJA_2022_pic<-subset(pic_sub,study=="McCullochJA_2022",drop=T)
SpencerCN_2021_pic<-subset(pic_sub,study=="SpencerCN_2021",drop=T)

LeeKA_2022_ratio<-dsv_pfs(LeeKA_2022_pic,"LeeKA_2022")
McCullochJA_2022_ratio<-dsv_pfs(McCullochJA_2022_pic,"McCullochJA_2022")
SpencerCN_2021_ratio<-dsv_pfs(SpencerCN_2021_pic,"SpencerCN_2021")


table(McCullochJA_2022_pic$temp)
p1<-ggplot(McCullochJA_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("McCullochJA_\n2022")+
   scale_fill_manual(breaks=c("PFS_less12","PFS_larger12"),labels=c("PFS<12m","PFS>=12m"),values = rev(mycolor2_yellow_green )) +
   labs(x = '', y = 'PFS_12_months (Ratio)')+
   geom_text(aes(y=new_col, label=value), size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=3)","Non-Delection (n=14)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=11,face="plain")
   )

table(SpencerCN_2021_pic$temp)
p2<-ggplot(SpencerCN_2021_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("SpencerCN_\n2021")+
   scale_fill_manual(breaks=c("PFS_less12","PFS_larger12"),labels=c("PFS<12m","PFS>=12m"),values = rev(mycolor2_yellow_green )) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=16)","Non-Delection (n=7)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size =12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

table(LeeKA_2022_pic$temp)
p3<-ggplot(LeeKA_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("LeeKA_\n2022")+
   scale_fill_manual(breaks=c("PFS_less12","PFS_larger12"),labels=c("PFS<12m","PFS>=12m"),values = rev(mycolor2_yellow_green )) +
   labs(x = '', y ='')+
   geom_text(aes(y=new_col, label=value),size=5,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=27)","Non-Delection(n=29)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size =12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_blank(),
   plot.title=element_text(size=11,face="plain")
   )

p_whole<-p1+p2+p3 +plot_layout(ncol=3) + plot_annotation(title = "Melenoma\ndSV Akkermansia muciniphila ATCC BAA-835:863~865 kbp\nMeta OR=0.22;Meta p=2.30e-3")

#p_whole<-plot_grid(p1,p2,p3,p4,
#              rel_widths = c(0.55,0.45,0.45,0.7),align = 'hv',
#              ncol = 4,label_size= 10,vjust = 0)

tiff(file = "pics/sites/F5_dsv_mel_pfs_A.muciniphila_863_865_barplot.tiff", width =1700, height =1300, res =300) 
p_whole
dev.off()


#############################################################
## pfs  Faecalibacterium cf. prausnitzii KLE1255:458_459

pic_dsv$temp<-pic_dsv$"Faecalibacterium cf. prausnitzii KLE1255:458_459"
lg_res_model <- glm(pfs_12_months~cohort+temp+age+read_count,family=binomial(logit), data=pic_dsv)
lg_res <- summary(lg_res_model)

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(pfs_12_months)==F,drop=T)
table(pic_sub$study)

LeeKA_2022_pic<-subset(pic_sub,study=="LeeKA_2022",drop=T)
McCullochJA_2022_pic<-subset(pic_sub,study=="McCullochJA_2022",drop=T)
SpencerCN_2021_pic<-subset(pic_sub,study=="SpencerCN_2021",drop=T)

LeeKA_2022_ratio<-dsv_pfs(LeeKA_2022_pic,"LeeKA_2022")
McCullochJA_2022_ratio<-dsv_pfs(McCullochJA_2022_pic,"McCullochJA_2022")
SpencerCN_2021_ratio<-dsv_pfs(SpencerCN_2021_pic,"SpencerCN_2021")

table(McCullochJA_2022_pic$temp)
p1<-ggplot(McCullochJA_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("McCullochJA_\n2022")+
   scale_fill_manual(breaks=c("PFS_less12","PFS_larger12"),labels=c("PFS<12m","PFS>=12m"),values = rev(mycolor2_yellow_green )) +
   labs(x = '', y = 'PFS_12_months (Ratio)')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=24)","Non-Delection(n=11)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=10,face="plain")
   )

table(SpencerCN_2021_pic$temp)
p2<-ggplot(SpencerCN_2021_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("SpencerCN_\n2021")+
   scale_fill_manual(breaks=c("PFS_less12","PFS_larger12"),labels=c("PFS<12m","PFS>=12m"),values = rev(mycolor2_yellow_green )) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=22)","Non-Delection(n=16)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=10,face="plain")
   )

table(LeeKA_2022_pic$temp)
p3<-ggplot(LeeKA_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("LeeKA_\n2022")+
   scale_fill_manual(breaks=c("PFS_less12","PFS_larger12"),labels=c("PFS<12m","PFS>=12m"),values = rev(mycolor2_yellow_green )) +
   labs(x = '', y ='')+
   geom_text(aes(y=new_col, label=value),size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=27)","Non-Delection(n=39)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_blank(),
   plot.title=element_text(size=10,face="plain")
   )

p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "Melenoma\ndSV Faecalibacterium cf. prausnitzii KLE1255:458 ~ 459 kbp")

#p_whole<-plot_grid(p1,p2,p3,p4,
#              rel_widths = c(0.55,0.45,0.45,0.7),align = 'hv',
#              ncol = 4,label_size= 10,vjust = 0)

tiff(file = "pics/sites/F5_dsv_mel_pfs_F.prausnitzii_458_459_barplot.tiff", width =1600, height =1500, res =300) 
p_whole
dev.off()


##################################################################################
### vsgv picture
all_basic$id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Coprococcus catus GD/7:648_649",
"Roseburia intestinalis XB6B4:2891_2892",
"Akkermansia muciniphila ATCC BAA-835:2546_2547",
"Akkermansia muciniphila ATCC BAA-835:2313_2314",
"Akkermansia muciniphila ATCC BAA-835:2314_2315",
"Akkermansia muciniphila ATCC BAA-835:2319_2320",
"Akkermansia muciniphila ATCC BAA-835:2546_2547",
"Akkermansia muciniphila ATCC BAA-835:2547_2548;2549_2550"
)]


vsgv_pic$id<-row.names(vsgv_pic)
### mel 
pic_vsv<-merge(mel_basic,vsgv_pic,by.x="id",by.y="id")

#################################################################################
## pfs Coprococcus catus GD/7:648_649

pic_vsv$temp<-pic_vsv$"Coprococcus catus GD/7:648_649"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(pfs_12_months)==F,drop=T)
table(pic_sub$study)

LeeKA_2022_pic<-subset(pic_sub,study=="LeeKA_2022",drop=T)
McCullochJA_2022_pic<-subset(pic_sub,study=="McCullochJA_2022",drop=T)

LeeKA_2022_ratio<-Q4_vsv_pfs(LeeKA_2022_pic,"LeeKA_2022\n(n=42)")
McCullochJA_2022_ratio<-Q4_vsv_pfs(McCullochJA_2022_pic,"McCullochJA_2022\n(n=11)")

ratio_pic<-rbind(McCullochJA_2022_ratio,LeeKA_2022_ratio)
ratio_pic<-data.frame(ratio_pic[seq(2,nrow(ratio_pic),2),])
ratio_pic$Qt<-as.numeric(ratio_pic$variable)

tiff(file = "pics/sites/F5_vsv_mel_pfs_Coprococcus_catus_648_649.tiff", width =1300, height =1100, res =300) 

ggplot(ratio_pic,aes(Qt,value,color=study))+
  geom_point(size=2)+
  geom_line(position = position_dodge(0.1),cex=1)+
  scale_x_continuous(breaks =seq(1,4),labels=c("Q1","Q2","Q3","Q4"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  geom_text(aes(label=value),size=4,vjust=-0.9)+
  ggtitle("Melenoma\nvSV Coprococcus catus GD/7:648~649 kbp\nMeta OR=0.38;Meta p=2.93e-3")+
  scale_color_brewer(palette = 'Set1')+
  labs(x='vSV',y='PFS >= 12 months (%)')+
  theme_pander()+
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(0.5, 'cm'),plot.title=element_text(size=12,face="plain"))

dev.off()



#################################################################################
## pfs Roseburia intestinalis XB6B4:2891_2892

pic_vsv$temp<-pic_vsv$"Roseburia intestinalis XB6B4:2891_2892"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(pfs_12_months)==F,drop=T)
table(pic_sub$study)

LeeKA_2022_pic<-subset(pic_sub,study=="LeeKA_2022",drop=T)
McCullochJA_2022_pic<-subset(pic_sub,study=="McCullochJA_2022",drop=T)
SpencerCN_2021_pic<-subset(pic_sub,study=="SpencerCN_2021",drop=T)
PetersBA_2019_pic<-subset(pic_sub,study=="PetersBA_2019",drop=T)

LeeKA_2022_ratio<-Q4_vsv_pfs(LeeKA_2022_pic,"LeeKA_2022\n(n=31)")
McCullochJA_2022_ratio<-Q4_vsv_pfs(McCullochJA_2022_pic,"McCullochJA_2022\n(n=23)")
#SpencerCN_2021_ratio<-Q4_vsv_pfs(SpencerCN_2021_pic,"SpencerCN_2021\n(n=4)")

ratio_pic<-rbind(McCullochJA_2022_ratio,LeeKA_2022_ratio)
ratio_pic<-data.frame(ratio_pic[seq(2,nrow(ratio_pic),2),])
ratio_pic$Qt<-as.numeric(ratio_pic$variable)

tiff(file = "pics/sites/F5_vsv_mel_pfs_Roseburia_intestinalis_1931_1932.tiff", width =1300, height =1300, res =300) 

ggplot(ratio_pic,aes(Qt,value,color=study))+
  geom_point(size=2)+
  geom_line(position = position_dodge(0.1),cex=1)+
  scale_x_continuous(breaks =seq(1,4),labels=c("Q1","Q2","Q3","Q4"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  geom_text(aes(label=value),size=4,vjust=-0.9)+
  ggtitle("Melenoma\nvSV R.intestinalis XB6B4:2891~2892 kbp\nMeta OR=0.45;Meta p=6.93e-3")+
  scale_color_brewer(palette = 'Set1')+
  labs(x='vSV',y='PFS >= 12 months (%)')+
  theme_pander()+
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(0.5, 'cm'),plot.title=element_text(size=12,face="plain"))

dev.off()



#################################################################################
## response Akkermansia muciniphila ATCC BAA-835:2546_2547

pic_vsv$temp<-pic_vsv$"Akkermansia muciniphila ATCC BAA-835:2546_2547"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(response_code)==F,drop=T)

table(pic_sub$study)

FrankelAE_2017_pic<-subset(pic_sub,study=="FrankelAE_2017",drop=T)
LeeKA_2022_pic<-subset(pic_sub,study=="LeeKA_2022",drop=T)
McCullochJA_2022_pic<-subset(pic_sub,study=="McCullochJA_2022",drop=T)
SpencerCN_2021_pic<-subset(pic_sub,study=="SpencerCN_2021",drop=T)

FrankelAE_2017_ratio<-Q4_vsv_resp(FrankelAE_2017_pic,"FrankelAE_2017\n(n=8)")
LeeKA_2022_ratio<-Q4_vsv_resp(LeeKA_2022_pic,"LeeKA_2022\n(n=56)")
McCullochJA_2022_ratio<-Q4_vsv_resp(McCullochJA_2022_pic,"McCullochJA_2022\n(n=17)")
SpencerCN_2021_ratio<-Q4_vsv_resp(SpencerCN_2021_pic,"SpencerCN_2021\n(n=30)")

ratio_pic<-rbind(SpencerCN_2021_ratio,LeeKA_2022_ratio,McCullochJA_2022_ratio)
ratio_pic<-data.frame(ratio_pic[seq(2,nrow(ratio_pic),2),])
ratio_pic$Qt<-as.numeric(ratio_pic$variable)

tiff(file = "pics/sites/F5_vsv_mel_response_Akkermansia muciniphila_2546_2547.tiff", width =1300, height =1300, res =300) 

ggplot(ratio_pic,aes(Qt,value,color=study))+
  geom_point(size=2)+
  geom_line(position = position_dodge(0.1),cex=1)+
  scale_x_continuous(breaks =seq(1,4),labels=c("Q1","Q2","Q3","Q4"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  geom_text(aes(label=value),size=4,vjust=-0.9)+
  ggtitle("Melenoma\nvSV A.muciniphila BAA-835:\n2546~2547 kbp\nmeta OR=1.86; meta p=5.04e-3")+
  scale_color_brewer(palette = 'Set1')+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  labs(x='vSV',y='Responser (CR/PR,%)')+
  theme_pander()+
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(0.5, 'cm'),plot.title=element_text(size=12,face="plain"))

dev.off()



#################################################################################
## response Akkermansia muciniphila ATCC BAA-835:2319_2320

pic_vsv$temp<-pic_vsv$"Akkermansia muciniphila ATCC BAA-835:2314_2315"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(response_code)==F,drop=T)

table(pic_sub$study)

FrankelAE_2017_pic<-subset(pic_sub,study=="FrankelAE_2017",drop=T)
LeeKA_2022_pic<-subset(pic_sub,study=="LeeKA_2022",drop=T)
McCullochJA_2022_pic<-subset(pic_sub,study=="McCullochJA_2022",drop=T)
SpencerCN_2021_pic<-subset(pic_sub,study=="SpencerCN_2021",drop=T)

FrankelAE_2017_ratio<-Q4_vsv_resp(FrankelAE_2017_pic,"FrankelAE_2017\n(n=8)")
LeeKA_2022_ratio<-Q4_vsv_resp(LeeKA_2022_pic,"LeeKA_2022\n(n=56)")
McCullochJA_2022_ratio<-Q4_vsv_resp(McCullochJA_2022_pic,"McCullochJA_2022\n(n=17)")
SpencerCN_2021_ratio<-Q4_vsv_resp(SpencerCN_2021_pic,"SpencerCN_2021\n(n=30)")

ratio_pic<-rbind(SpencerCN_2021_ratio,LeeKA_2022_ratio,McCullochJA_2022_ratio)
ratio_pic<-data.frame(ratio_pic[seq(2,nrow(ratio_pic),2),])
ratio_pic$Qt<-as.numeric(ratio_pic$variable)

tiff(file = "pics/sites/F5_vsv_mel_response_Akkermansia muciniphila_2514_2515.tiff", width =1350, height =1100, res =300) 

ggplot(ratio_pic,aes(Qt,value,color=study))+
  geom_point(size=2)+
  geom_line(position = position_dodge(0.1),cex=1)+
  scale_x_continuous(breaks =seq(1,4),labels=c("Q1","Q2","Q3","Q4"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  geom_text(aes(label=value),size=4,vjust=-0.9)+
  ggtitle("Melenoma\nvSV A.muciniphila BAA-835:\n2514 ~ 2515 kbp\nmeta OR=1.98; meta p=1.49e-3")+
  scale_color_brewer(palette = 'Set1')+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  labs(x='vSV',y='Responser (CR/PR,%)')+
  theme_pander()+
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(0.5, 'cm'),plot.title=element_text(size=12,face="plain"))

dev.off()


