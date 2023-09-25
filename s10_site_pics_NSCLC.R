### Survival and barplot for sites
### 2023-09-21
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
library("patchwork")

source("functions.R")

###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/France_Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/France_vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/France_dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")

NSCLC_basic<-subset(all_basic,cancer_type=="NSCLC",drop=T)
RCC_basic<-subset(all_basic,cancer_type=="RCC",drop=T)

load("01.cleanData/SV_all/France_dsgv.RData")
load("01.cleanData/SV_all/France_vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub

#######################################################################################
### clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Akkermansia muciniphila ATCC BAA-835:731_734",
"Alistipes obesi:1922_1923",
"Bacteroides caccae ATCC 43185:0_1;1626_1628",
"Bacteroides caccae ATCC 43185:1713_1715",
"Bacteroides caccae ATCC 43185:1990_1991",
"Bacteroides xylanisolvens XB1A:638_641",
"Faecalibacterium cf. prausnitzii KLE1255:1563_1564",
"Parabacteroides distasonis ATCC 8503:4655_4656",
"Prevotella copri DSM 18205:753_754;2664_2667",
"Ruminococcus bicirculans:1181_1182",
"Ruminococcus bicirculans:284_286",
"Ruminococcus bromii L2-63:1083_1084",
"Akkermansia muciniphila ATCC BAA-835:1624_1626",
"Akkermansia muciniphila ATCC BAA-835:2528_2529;2532_2550",
"Akkermansia muciniphila ATCC BAA-835:83_84",
"Coprococcus comes ATCC 27758:1400_1402"
)]


########################################################################################
##  response (NSCLC)

dsgv_pic$id<-row.names(dsgv_pic)
pic_dsv<-merge(NSCLC_basic,dsgv_pic,by.x="id",by.y="id")

########################################################################################
##  response Akkermansia muciniphila ATCC BAA-835:731_734

pic_dsv$temp<-pic_dsv$"Akkermansia muciniphila ATCC BAA-835:731_734"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(response_code)==F,drop=T)
table(pic_sub$study)

DerosaL_2022_pic<-subset(pic_sub,study=="DerosaL_2022",drop=T)
RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)

DerosaL_2022_ratio<-dsv_resp(DerosaL_2022_pic,"DerosaL_2022")
RoutyB_2018_ratio<-dsv_resp(RoutyB_2018_pic,"RoutyB_2018")

table(DerosaL_2022_pic$temp)

p1<-ggplot(DerosaL_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9)+ggtitle("DerosaL_\n2022")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=38)","Non-Delection (n=58)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size =12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=11,face="plain"),axis.text.y = element_text(size = 12),
   )

table(RoutyB_2018_pic$temp)

p2<-ggplot(RoutyB_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9)+ggtitle("RoutyB_\n2018")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection (n=15)","Non-Delection (n=26)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size =12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),axis.text.y = element_blank(),
   plot.title=element_text(size=11,face="plain")
   )

p_whole<-p1+p2+plot_layout(ncol=2) + plot_annotation(title = "NSCLC\ndSV A.muciniphila BAA-835:731~734 kbp\nMeta OR=4.50; Meta p=2.30e-3")

tiff(file = "pics/sites/F6_dsv_resp_NSCLC_Akk_731_734_barplot.tiff", width =1200, height =1500, res =300) 
p_whole
dev.off()



########################################################################################
##  response Parabacteroides distasonis ATCC 8503:4655_4656

pic_dsv$temp<-pic_dsv$"Parabacteroides distasonis ATCC 8503:4655_4656"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(response_code)==F,drop=T)
table(pic_sub$study)

DerosaL_2022_pic<-subset(pic_sub,study=="DerosaL_2022",drop=T)
RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)

DerosaL_2022_ratio<-dsv_resp(DerosaL_2022_pic,"DerosaL_2022")
RoutyB_2018_ratio<-dsv_resp(RoutyB_2018_pic,"RoutyB_2018")

table(DerosaL_2022_pic$temp)

p1<-ggplot(DerosaL_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9)+ggtitle("DerosaL_\n2022")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=68)","Non-Delection(n=172)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=11,face="plain"),axis.text.y = element_text(size = 12),
   )

table(RoutyB_2018_pic$temp)

p2<-ggplot(RoutyB_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9)+ggtitle("RoutyB_\n2018")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=22)","Non-Delection(n=69)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12),axis.text.y = element_blank(),
   plot.title=element_text(size=11,face="plain")
   )

p_whole<-p1+p2 +plot_layout(ncol=2) + plot_annotation(title = "NSCLC\ndSV P.distasonis ATCC 8503:4655~4656 kbp\nMeta OR=0.39; Meta p=8.94e-4")

tiff(file = "pics/sites/F6_dsv_resp_NSCLC_Pdistasonis_4655_4656_barplot.tiff", width =1200, height =1500, res =300) 
p_whole
dev.off()




#################################################################################
## os Akkermansia muciniphila ATCC BAA-835:2528_2529;2532_2550
pic_sub<-pic_dsv
pic_sub$subtype<-pic_dsv$"Akkermansia muciniphila ATCC BAA-835:2528_2529;2532_2550"
a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_NSCLC_OS_Akkermansiamuciniphila2528_2529_2532_2550.tiff", width =1300, height =1500, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = mycolor2_green_blue, 
conf.int = F, pval ="log rank P=0.002",  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(0,0.15),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="NSCLC\nA.muciniphila BAA-835:\n2528~2529;2532~2550 kbp",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Delection","Non-delection"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()



#################################################################################
## os Coprococcus comes ATCC 27758:1400_1402
pic_sub<-pic_dsv
pic_sub$subtype<-pic_sub$"Coprococcus comes ATCC 27758:1400_1402"
a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_NSCLC_OS_Coprococcuscomes_1400_1402.tiff", width =1300, height =1500, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = mycolor2_green_blue, 
conf.int = F, pval ="log rank P=1.8e-4",,  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(0,0.1),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="NSCLC\nC.comes ATCC 27758:1400~1402 kbp",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Non-delection","Delection"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()



##################################################################################
### vsgv picture
all_basic$id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Bacteroides caccae ATCC 43185:148_150 and 8 segments",
"Bacteroides caccae ATCC 43185:1645_1646 and 10 segments",
"Bacteroides caccae ATCC 43185:1992_1998",
"Bacteroides caccae ATCC 43185:3235_3238",
"Blautia wexlerae DSM 19850:659_660",
"Ruminococcus lactaris ATCC 29176:570_571",
"Bacteroides clarus YIT 12056:3156_3159",
"Ruminococcus bromii L2-63:2142_2143 and 5 segments",
"Ruminococcus bromii L2-63:2192_2194 and 2 segments",
"Ruminococcus bromii L2-63:40_45 and 8 segments"
)]


vsgv_pic$id<-row.names(vsgv_pic)

### NSCLC
pic_vsv<-merge(NSCLC_basic,vsgv_pic,by.x="id",by.y="id")

#################################################################################
## os Ruminococcus bromii L2-63:40_45 and 8 segments

library(RColorBrewer)
#green_colors<-brewer.pal(6,"Greens")
green_colors<-c("#A1D99B","#74C476","#31A354","#006D2C")

pic_sub<-pic_vsv
pic_sub$temp<-pic_sub$"Ruminococcus bromii L2-63:40_45 and 8 segments"

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_NSCLC_os_vsv_Rbromii_40_45.tiff", width =1600, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = green_colors, 
conf.int = F, pval ="log rank P=0.0055", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(0,0.15),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="NSCLC:\nR.bromii L2-63:40~45 and 8 segments kbp",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Q1","Q2","Q3","Q4"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()


#################################################################################
## os Bacteroides clarus YIT 12056:3156_3159

library(RColorBrewer)
#green_colors<-brewer.pal(6,"Greens")
green_colors<-c("#A1D99B","#74C476","#31A354","#006D2C")

pic_sub<-pic_vsv
pic_sub$temp<-pic_sub$"Bacteroides clarus YIT 12056:3156_3159"

qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
 pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }

a<-summary(coxph(Surv(os, os_event)~subtype,pic_sub,na.action=na.exclude))

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~subtype, data =pic_sub) 
print(summary(coxph(Surv(os, os_event) ~subtype, pic_sub)))

tiff(file = "pics/sites/F6_NSCLC_os_vsv_Bacteroides_clarus_3156_3159.tiff", width =1600, height =1600, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 1, 
palette = green_colors, 
conf.int = F, pval ="log rank P=0.0013",  
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(0,0.15),font.title=14,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="NSCLC:\nB.clarus YIT 12056:3156~3159 kbp",
xlab="Time(year)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Q1","Q2","Q3","Q4"),
legend.title="Subtype",risk.table.height = 0.35,surv.plot.height = 0.65
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()



### vsv NSCLC response 

#################################################################################
## response 
## Bacteroides caccae ATCC 43185:148_150 and 8 segments
## Ruminococcus lactaris ATCC 29176:2420_2421
## Ruminococcus lactaris ATCC 29176:570_571

pic_vsv$temp<-pic_vsv$"Bacteroides caccae ATCC 43185:148_150 and 8 segments"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(response_code)==F,drop=T)

table(pic_sub$study)

DerosaL_2022_pic<-subset(pic_sub,study=="DerosaL_2022",drop=T)
RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)

DerosaL_2022_ratio<-Q4_vsv_resp(DerosaL_2022_pic,"DerosaL_2022\n(n=148)")
RoutyB_2018_ratio<-Q4_vsv_resp(RoutyB_2018_pic,"RoutyB_2018\n(n=60)")

ratio_pic<-rbind(DerosaL_2022_ratio,RoutyB_2018_ratio)
ratio_pic<-data.frame(ratio_pic[seq(2,nrow(ratio_pic),2),])
ratio_pic$Qt<-as.numeric(ratio_pic$variable)

tiff(file = "pics/sites/F6_vsv_NSCLC_response_Bcaccae_148_150.tiff", width =1200, height =1100, res =300) 

ggplot(ratio_pic,aes(Qt,value,color=study))+
  geom_point(size=2)+
  geom_line(position = position_dodge(0.1),cex=1)+
  scale_x_continuous(breaks =seq(1,4),labels=c("Q1","Q2","Q3","Q4"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  ggtitle("NSCLC\nvSV B.caccae ATCC 43185:148_150 and\n8 segments kbp\nMeta OR=1.49; Meta p=8.63e-3")+
  geom_text(aes(label=value),size=4,vjust=-0.9)+
  scale_color_brewer(palette = 'Set1')+
  labs(x='vSV',y='Responser (%)')+
  theme_pander()+
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(0.5, 'cm'),plot.title=element_text(size=12,face="plain"))

dev.off()



#################################################################################
## response 
## Ruminococcus lactaris ATCC 29176:570_571

pic_vsv$temp<-pic_vsv$"Ruminococcus lactaris ATCC 29176:570_571"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(response_code)==F,drop=T)

table(pic_sub$study)

DerosaL_2022_pic<-subset(pic_sub,study=="DerosaL_2022",drop=T)
RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)

DerosaL_2022_ratio<-Q4_vsv_resp(DerosaL_2022_pic,"DerosaL_2022\n(n=61)")
RoutyB_2018_ratio<-Q4_vsv_resp(RoutyB_2018_pic,"RoutyB_2018\n(n=22)")

ratio_pic<-rbind(DerosaL_2022_ratio,RoutyB_2018_ratio)
ratio_pic<-data.frame(ratio_pic[seq(2,nrow(ratio_pic),2),])
ratio_pic$Qt<-as.numeric(ratio_pic$variable)

tiff(file = "pics/sites/F6_vsv_NSCLC_response_RlactarisATCC29176570_571.tiff", width =1200, height =1100, res =300) 

ggplot(ratio_pic,aes(Qt,value,color=study))+
  geom_point(size=2)+
  geom_line(position = position_dodge(0.1),cex=1)+
  scale_x_continuous(breaks =seq(1,4),labels=c("Q1","Q2","Q3","Q4"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  ggtitle("NSCLC\nvSV R.lactaris ATCC 29176:570~571 kbp\n meta OR=2.70; meta P=8.1e-4")+
  geom_text(aes(label=value),size=4,vjust=-0.9)+
  scale_color_brewer(palette = 'Set1')+
  labs(x='vSV',y='Responser (%)')+
  theme_pander() +
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(0.5, 'cm'),plot.title=element_text(size=12,face="plain"))

dev.off()


#################################################################################
## response 
## Blautia wexlerae DSM 19850:659_660

pic_vsv$temp<-pic_vsv$"Blautia wexlerae DSM 19850:659_660"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(response_code)==F,drop=T)

table(pic_sub$study)

DerosaL_2022_pic<-subset(pic_sub,study=="DerosaL_2022",drop=T)
RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)

DerosaL_2022_ratio<-Q4_vsv_resp(DerosaL_2022_pic,"DerosaL_2022\n(n=92)")
RoutyB_2018_ratio<-Q4_vsv_resp(RoutyB_2018_pic,"RoutyB_2018\n(n=20)")

ratio_pic<-rbind(DerosaL_2022_ratio,RoutyB_2018_ratio)
ratio_pic<-data.frame(ratio_pic[seq(2,nrow(ratio_pic),2),])
ratio_pic$Qt<-as.numeric(ratio_pic$variable)

tiff(file = "pics/sites/F6_vsv_NSCLC_response_Blautia_wexlerae_659_660.tiff", width =1200, height =1000, res =300) 

ggplot(ratio_pic,aes(Qt,value,color=study))+
  geom_point(size=2)+
  geom_line(position = position_dodge(0.1),cex=1)+
  scale_x_continuous(breaks =seq(1,4),labels=c("Q1","Q2","Q3","Q4"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  ggtitle("NSCLC\nvSV B.wexlerae DSM 19850:659_660")+
  geom_text(aes(label=value),size=4,vjust=-0.9)+
  scale_color_brewer(palette = 'Set1')+
  labs(x='vSV',y='Responser (%)')+
  theme_pander() +
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(0.5, 'cm'),plot.title=element_text(size=12,face="plain"))

dev.off()



############################################################################################
### RCC

#######################################################################################
### clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Alistipes putredinis DSM 17216:107_109",
"Alistipes putredinis DSM 17216:118_120 and 2 segments",
"Alistipes putredinis DSM 17216:123_124",
"Parabacteroides distasonis ATCC 8503:2192_2195"
)]


########################################################################################
##  response (RCC)

dsgv_pic$id<-row.names(dsgv_pic)
pic_dsv<-merge(RCC_basic,dsgv_pic,by.x="id",by.y="id")

########################################################################################
##  response "Alistipes putredinis DSM 17216:107_109"

pic_dsv$temp<-pic_dsv$"Alistipes putredinis DSM 17216:107_109"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(response_code)==F,drop=T)
table(pic_sub$study)

RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)
RoutyB_2018_ratio<-dsv_resp(RoutyB_2018_pic,"RoutyB_2018")

table(RoutyB_2018_pic$temp)

p1<-ggplot(RoutyB_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9)+ggtitle("RoutyB_2018")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=48)","Non-Delection(n=24)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size =12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans"),legend.position="right",
   plot.title=element_text(size=11,face="plain"),axis.text.y = element_text(size = 12)
   )

p_whole<-p1+plot_layout(ncol=1) + plot_annotation(title = "RCC\ndSV A.putredinis DSM 17216:107~109 kbp \n Meta OR=7.44; Meta p=4.03e-3")

tiff(file = "pics/sites/F7_dsv_resp_RCC_A.putredinis_107_109_barplot.tiff", width =900, height =1500, res =300) 
p_whole
dev.off()


########################################################################################
##  response "Alistipes putredinis DSM 17216:118_120 and 2 segments"

pic_dsv$temp<-pic_dsv$"Alistipes putredinis DSM 17216:118_120 and 2 segments"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(response_code)==F,drop=T)
table(pic_sub$study)

RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)
RoutyB_2018_ratio<-dsv_resp(RoutyB_2018_pic,"RoutyB_2018")

table(RoutyB_2018_pic$temp)

p1<-ggplot(RoutyB_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9)+ggtitle("RoutyB_2018")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=48)","Non-Delection(n=24)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans"),legend.position="right",
   plot.title=element_text(size=12,face="plain")
   )

p_whole<-p1+plot_layout(ncol=1) + plot_annotation(title = "NSCLC\ndSV A.putredinis DSM 17216:118_120 and 2 segments")

tiff(file = "pics/sites/F7_dsv_resp_RCC_A.putredinis_118_120and2segments_barplot.tiff", width =900, height =1500, res =300) 
p_whole
dev.off()



########################################################################################
##  response Alistipes putredinis DSM 17216:123_124

pic_dsv$temp<-pic_dsv$"Alistipes putredinis DSM 17216:123_124"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(response_code)==F,drop=T)
table(pic_sub$study)

RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)
RoutyB_2018_ratio<-dsv_resp(RoutyB_2018_pic,"RoutyB_2018")

table(RoutyB_2018_pic$temp)

p1<-ggplot(RoutyB_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9)+ggtitle("RoutyB_2018")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=49)","Non-Delection(n=23)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 12)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans"),legend.position="right",
   plot.title=element_text(size=12,face="plain"),axis.text.y = element_text(size = 12)
   )

p_whole<-p1+plot_layout(ncol=1) + plot_annotation(title = "NSCLC\ndSV A.putredinis DSM 17216:123~124 kbp")

tiff(file = "pics/sites/F7_dsv_resp_RCC_A.putredinis_123_124_barplot.tiff", width =800, height =1300, res =300) 
p_whole
dev.off()


###########################################################################################################
## response Parabacteroides distasonis ATCC 8503:2192_2195

pic_dsv$temp<-pic_dsv$"Parabacteroides distasonis ATCC 8503:2192_2195"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(response_code)==F,drop=T)
table(pic_sub$study)

RoutyB_2018_pic<-subset(pic_sub,study=="RoutyB_2018",drop=T)
RoutyB_2018_ratio<-dsv_resp(RoutyB_2018_pic,"RoutyB_2018")

table(RoutyB_2018_pic$temp)

p1<-ggplot(RoutyB_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9)+ggtitle("RoutyB_2018")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=4,colour="White")+
   scale_x_discrete(limits=c("Delection","Non-delection"),
   labels=c("Delection(n=64)","Non-Delection(n=18)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 8)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 45, hjust = 1),text=element_text(family="sans"),legend.position="right",
   plot.title=element_text(size=8,face="bold")
   )

p_whole<-p1+plot_layout(ncol=1) + plot_annotation(title = "NSCLC\ndSV P.distasonis:2192_2195")

tiff(file = "pics/sites/F7_dsv_resp_RCC_P.distasonis_2192_2195_barplot.tiff", width =900, height =1500, res =300) 
p_whole
dev.off()




##### vsv



"Barnesiella intestinihominis YIT 11860:333_335"




