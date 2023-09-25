###################################################################################################
####### draw pictures (Figure 1)
####### draw starked barplot for response  
####### Liurong
####### 2023-09-18

##install.packages("car")

library("ggplot2")
library("reshape2")	
library("RColorBrewer")
#library("plyr")

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")
source("functions.R")

# Basic
France_basic<-read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")
USA_UK_basic<-read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv")
clin<-rbind(France_basic,USA_UK_basic)
write.table(clin,"01.cleanData/phen/Clinical_basic_overlap.tsv",sep = "\t", quote = F, col.names = T, row.names = F)

table(clin$cancer_type)
table(clin$gender)

table(clin$study,clin$target)
table(clin$study,clin$cancer_type)

# clin<-read.csv("00.rawData/clinical/clinical_whole.csv",header=T)
# table(clin$cancer_type)

### basic characteristic 
pfs_basic<-subset(clin,pfs!="NA" & pfs_event!="NA",drop=T)
nrow(pfs_basic)
table(pfs_basic$study)
table(pfs_basic$cancer_type)

pfs_6m<-subset(clin,pfs_6_months!="NA",drop=T)
nrow(pfs_6m)

table(pfs_6m$cohort)
table(pfs_6m$cancer_type)
table(pfs_6m$study)

table(clin$pfs_6_months)
table(clin$pfs_12_months)


##########################################################################
##  clinical characteristic 
## Gender

all_gender_tbl <- table(clin$study, clin$gender)
chisq.test(all_gender_tbl) 
clin$Cohort<-as.factor(clin$study)

nr<-nrow(clin)
for(i in 1:nr){
  if(is.na(clin$gender[i])==T){clin$Gender[i]="Unknown"}
  if(is.na(clin$gender[i])==F){
    if (clin$gender[i]=="female"){clin$Gender[i]="female"}
    else if (clin$gender[i]=="male"){clin$Gender[i]="male"}
  }
 }

clin$Gender<-as.factor(clin$Gender)

##par(family = "sans")

p_gender<-ggplot(data=clin)+
  geom_mosaic(aes(x = product(study), fill=Gender))+
  ylab("Gender")+
  xlab("Study")+
  scale_fill_manual(values=mycolor2_yellow_green) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.ticks.length = unit(0, "cm"), 
        axis.text.x = element_text(angle = 60, hjust = 1,colour = "black",size=15),
        axis.text.y = element_text(colour = "black",size=15), 
        axis.title = element_text(colour = "black",size=18), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_gender.tiff", width =2000, height =1800, res =300) 
p_gender
dev.off()



## Age
mean(as.numeric(clin$age), na.rm = T)
se(as.numeric(clin$age))
min(as.numeric(clin$age), na.rm = T)
max(as.numeric(clin$age), na.rm = T)

clin$age<-as.numeric(clin$age)

## annova test 
##wilcox.test(clin$age~clin$study)

clin_age<-subset(clin,age!="NA",drop=T)

p_age_cont<-ggplot(clin_age,aes(x=age, color = dataset, fill =dataset))+
  geom_density(alpha = 0.2)+
  geom_rug(alpha = 0.5,length = unit(0.05, "npc"))+
  ylab('Density')+
  xlab('Age (year)')+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(name=NULL,
                    breaks=c("PRJEB22863","PRJNA397906","PRJNA541981","PRJNA762360","PRJNA770295","PRJNA751792"),
                    labels=c("RoutyB_2018","FrankelAE_2017","PetersBA_2019","McCullochJA_2022","SpencerCN_2021","DerosaL_2022"),
                     values = mycolor6)+
  scale_fill_manual(name=NULL,
                    breaks=c("PRJEB22863","PRJNA397906","PRJNA541981","PRJNA762360","PRJNA770295","PRJNA751792"),
                    labels=c("RoutyB_2018","FrankelAE_2017","PetersBA_2019","McCullochJA_2022","SpencerCN_2021","DerosaL_2022"),
                    values = mycolor6)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "bottom",
        legend.text=element_text(colour = "black",size=14),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=20),
        axis.text.y = element_text(colour = "black",size=20), 
        axis.title = element_text(colour = "black",size=18), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_age_num.tiff", width =2600, height =2000, res =300) 
p_age_cont
dev.off()



## Age

## annova test 
##wilcox.test(clin$age~clin$study)
clin_age<-subset(clin,age_bins!="NA",drop=T)

p_age<-ggplot(clin_age,aes(x=age_bins, color = dataset, fill =dataset))+
  geom_density(alpha = 0.2)+
  geom_rug(alpha = 0.5,length = unit(0.05, "npc"))+
  #geom_histogram(aes(y=after_stat(count / sum(count)),
  #                   fill=dataset),stat="count",alpha=0.5,color="black")+
  ylab('Density')+
  xlab('Age (year)')+
  scale_x_continuous(breaks = c(4,8,12,16),labels=c("35-40","55-60","70-75","90-95"))+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(name=NULL,
                    breaks=c("PRJEB22863","PRJNA397906","PRJNA541981","PRJNA762360","PRJNA770295","PRJNA751792","PRJEB43119"),
                    labels=c("RoutyB_2018","FrankelAE_2017","PetersBA_2019","McCullochJA_2022","SpencerCN_2021","DerosaL_2022","LeeKA_2022"),
                     values = mycolor7)+
  scale_fill_manual(name=NULL,
                    breaks=c("PRJEB22863","PRJNA397906","PRJNA541981","PRJNA762360","PRJNA770295","PRJNA751792","PRJEB43119"),
                    labels=c("RoutyB_2018","FrankelAE_2017","PetersBA_2019","McCullochJA_2022","SpencerCN_2021","DerosaL_2022","LeeKA_2022"),
                    values = mycolor7)+
  guides(col = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "bottom",
        legend.text=element_text(colour = "black",size=16),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_age_bins.tiff", width =2300, height =1800, res =300) 
p_age
dev.off()

########################################################################################
#### BMI
mean(as.numeric(clin$BMI), na.rm = T)
se(as.numeric(clin$BMI))
min(as.numeric(clin$BMI), na.rm = T)
max(as.numeric(clin$BMI), na.rm = T)

clin$BMI<-as.numeric(clin$BMI)
##wilcox.test(clin$BMI~clin$study)

clin_BMI<-subset(clin,BMI!="NA",drop=T)
table(clin_BMI$study)

p_bmi<-ggplot(clin_BMI,aes(x=BMI, color = dataset, fill =dataset))+
  geom_density(alpha = 0.2)+
  geom_rug(alpha = 0.5,length = unit(0.05, "npc"))+
  ylab('Density')+
  xlab('Body mass index')+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(name=NULL,
                    breaks=c("PRJNA762360","PRJNA770295","PRJEB43119","PRJNA541981"),
                    labels=c("McCullochJA_2022","SpencerCN_2021","LeeKA_2022","PetersBA_2019"),
                     values = c("#20854EFF","#7876B1FF","#FFDC91FF","#E18727FF"))+
  scale_fill_manual(name=NULL,
                    breaks=c("PRJNA762360","PRJNA770295","PRJEB43119","PRJNA541981"),
                    labels=c("McCullochJA_2022","SpencerCN_2021","LeeKA_2022","PetersBA_2019"),
                    values = c("#20854EFF","#7876B1FF","#FFDC91FF","#E18727FF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "bottom",
        legend.text=element_text(colour = "black",size=16),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_BMI.tiff",width =2200, height =1700, res =300) 
p_bmi
dev.off()



#mycolor6
#hist(mtcars$mpg,
#     breaks = 15,
#     col ="#20854EFF",
#)
#"#6F99ADFF","#E18727FF","#BC3C29FF","#7876B1FF","#0072B5FF",#20854EFF



library(cowplot)
plot_grid(p_gender,p_age, p_bmi,
          labels = LETTERS[1:6],
          ncol = 3) 

library(patchwork)

p_gender+p_age+p_bmi+plot_layout(ncol = 3)+plot_annotation(tag_levels = LETTERS[1:6]) 


## plot
p_title <- ggdraw() + 
    draw_label(
      'Cohort characteristics',
      fontface = 'bold', x = 0, hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))

p_demo_grid<-plot_grid(
    p_title, 
    plot_grid(p_gender,p_age, p_bmi,
              rel_widths = c(0.2, 1, 1),align = 'hv',
              labels = c("A", "B", "C"),
              ncol = 1,label_size	= 8,vjust = 0),
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1.5)
  )
  
  
print(p_demo_grid)

##if(!dir.exists("02.summary_statistics")){dir.create("02.summary_statistics")}
##pdf("pics/Figure_characteristic_basic_phen_diff.pdf", height = 6, width = 3)

tiff(file = "pics/F1_test_characteristic_basic.tiff", width =1500, height =2000, res =300) 
print(p_demo_grid)
dev.off()



##################################################################################
### response
clin<- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv",header=T)

clin_resp<-subset(clin,response_code!="NA",drop=T)

table(clin_resp$cancer_type,clin_resp$target)
clin_resp$cancer_target<-paste(clin_resp$cancer_type,clin_resp$target)
table(clin_resp$cancer_target)

mel<-subset(clin_resp,cancer_type=="melanoma",drop=T)
ccr<-subset(clin_resp,cancer_type=="RCC",drop=T)
NSCLC<-subset(clin_resp,cancer_type=="NSCLC",drop=T)

table(mel$response_code)/sum(table(mel$response_code))
table(ccr$response_code)/sum(table(ccr$response_code))
table(NSCLC$response_code)/sum(table(NSCLC$response_code))


############################################################################
## ratio plot (use R 3.6.2) 
ratio<-table(clin_resp$response_code,clin_resp$cancer_type)
table(clin_resp$cancer_type)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,3)

for(i in 1:3){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("Melanoma","NSCLC","Cell renal cell carcinoma")
rownames(ratio_new)<-c("Non_RE","RE")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

tiff(file = "pics/F1_clinical_response.tiff", width =1800, height =2600, res =300) 

ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("")+
   scale_fill_manual(breaks=c("Non_RE","RE"),labels=c("SD/PD","CR/PR"),values = rev(mycolor2_blue_red)) +
   labs(x = '', y = 'Response rate')+
   geom_text(aes(y=new_col, label=value), size=8,colour="White")+
   scale_x_discrete(limits=flevels,
   labels=c("Melanoma (n=415)","NSCLC (n=447)","RCC (n=101)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 22,colour="black"),
   axis.text = element_text(size = 22,colour="black"), 
   axis.title = element_text(size = 22,colour="black"), 
   legend.title = element_blank(), 
   legend.text = element_text(size = 22,colour="black"),plot.margin = unit(c(0.5,0.5,1,1), "cm"),
   legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1,colour="black"))

dev.off()



##################################################################################
### pfs longer than 12 months 

clin_pfs<-subset(clin,pfs_12_months!="NA",drop=T)
table(clin_pfs$cancer_type)
table(clin_pfs$study)

table(mel$pfs_12_months)/sum(table(mel$pfs_12_months))

LeeKA_2022<-subset(clin_pfs,study=="LeeKA_2022",drop=T)
McCullochJA_2022<-subset(clin_pfs,study=="McCullochJA_2022",drop=T)
PetersBA_2020<-subset(clin_pfs,study=="PetersBA_2020",drop=T)
SpencerCN_2021<-subset(clin_pfs,study=="SpencerCN_2021",drop=T)

table(LeeKA_2022$pfs_12_months)/sum(table(LeeKA_2022$pfs_12_months))
table(McCullochJA_2022$pfs_12_months)/sum(table(McCullochJA_2022$pfs_12_months))
table(PetersBA_2020$pfs_12_months)/sum(table(PetersBA_2020$pfs_12_months))
table(SpencerCN_2021$pfs_12_months)/sum(table(SpencerCN_2021$pfs_12_months))


############################################################################
## ratio plot 
ratio<-table(clin_pfs$pfs_12_months,clin_pfs$study)
table(clin_pfs$study)

sum_vec<-apply(ratio,2,sum)
ratio_new<-matrix(NA,2,4)

for(i in 1:4){
 ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
 }

colnames(ratio_new)<-c("LeeKA_2022","McCullochJA_2022","SpencerCN_2021","PetersBA_2019")
rownames(ratio_new)<-c("pfs_less","pfs_longer")

#调整数据布局
ratio_pic<-data.frame(t(ratio_new))
ratio_pic$chara<-rownames(ratio_pic)

ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
ratio_pic<-melt(ratio_pic,id="chara")

flevels <- levels(ratio_pic$chara)
flevels <- rev(flevels)

ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic

tiff(file = "pics/F1_clinical_PFS.tiff", width =2500, height =2600, res =300) 

ggplot(ratio_pic, aes(chara, as.numeric(value), fill = variable)) +
  geom_col(position = 'stack', width =0.9) +ggtitle("")+
  scale_fill_manual(breaks=c("pfs_less","pfs_longer"),labels=c("PFS < 12 months","PFS ≥ 12 months"),values =c("#EE4C97FF","#0072B5FF")) +
  labs(x = '', y = 'PFS12 (ratio)')+
  geom_text(aes(y=new_col, label=value), size=10,colour="White")+
  scale_x_discrete(limits=flevels,
  labels=c("LeeKA_2022 (n=162)","McCullochJA_2022 (n=62)","SpencerCN_2021 (n=107)","PetersBA_2019 (n=21)"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
  strip.text = element_text(size = 22,colour="black"),
  axis.text = element_text(size = 22,colour="black"), 
  axis.title = element_text(size = 22,colour="black"), 
  legend.title = element_blank(), 
  legend.text = element_text(size = 22,colour="black"),plot.margin = unit(c(0.5,0.5,1,1), "cm"),
  legend.position="bottom",axis.text.x = element_text(angle = 45, hjust = 1,colour="black"))

dev.off()


######################################################################
### Survival plot between response and non-response group

library("survival")
library("ggplot2")
##install.packages('xfun')
library("survminer")

##########################################################################
### for the p values and survival
##os_1m_cancer<-subset(os_1m,cancer_type=="Melanoma" | cancer_type=="Urothelial cancer" | cancer_type=="Clear cell renal cell carcinoma",drop=T)

a<-summary(coxph(Surv(os, os_event)~response_code,clin,na.action=na.exclude))
table(clin$response_code)

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~response_code, data =clin) 
print(summary(coxph(Surv(os, os_event) ~response_code, clin)))

tiff(file = "pics/F1_OS_response.tiff", width =2500, height =3000, res =600) 

ggsurv<-ggsurvplot(yg.surv, size = 0.4, 
palette = rev(mycolor2_blue_red), 
conf.int = F, pval ="log rank p<2e-16",  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(2,0.9),font.title=24,
break.time.by =1,risk.table.fontsize=4,font.legend=14,
title="",
xlab="Time (years)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("SD/PD","CR/PR"),
legend.title="Subtype",risk.table.height = 0.3,surv.plot.height = 0.7)

ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()


a<-summary(coxph(Surv(os, os_event)~cancer_type,clin,na.action=na.exclude))

clin_os<-subset(clin,os!="NA",drop=T)
table(clin_os$cancer_type)

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(os, os_event) ~cancer_type, data =clin) 
print(summary(coxph(Surv(os, os_event) ~cancer_type, clin)))

tiff(file = "pics/F1_OS_cancer_type.tiff", width =2500, height =3000, res =600) 

ggsurv<-ggsurvplot(yg.surv, size = 0.4, 
palette =c("#E18727FF","#7876B1FF"), 
conf.int = F, pval ="log rank p=3e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(2,0.9),font.title=14,
break.time.by =1,risk.table.fontsize=4,font.legend=14,
title="",
xlab="Time (years)",ylab="Overall survival,%",risk.table.title="",
legend.labs=c("Melanoma","NSCLC"),
legend.title="Subtype",risk.table.height = 0.3,surv.plot.height = 0.7)

ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()


a<-summary(coxph(Surv(pfs, pfs_event)~response_code,clin,na.action=na.exclude))
table(clin$response_code)

###############################################################################################
### beautiful survival curves
yg.surv <- survfit(Surv(pfs, pfs_event) ~response_code, data =clin) 
print(summary(coxph(Surv(pfs, pfs_event) ~response_code, clin)))

tiff(file = "pics/FS1_PFS_response.tiff", width =1500, height =1500, res =300) 

ggsurv<-ggsurvplot(yg.surv, size = 0.4, 
palette = rev(mycolor2_blue_red), 
conf.int = F, pval ="log rank p<2e-16",  #"p=1e-6", 
risk.table = TRUE, 
risk.table.col = "strata", legend="none",pval.coord=c(2,0.9),font.title=16,
break.time.by =1,risk.table.fontsize=5,font.legend=12,
title="",
xlab="Time (years)",ylab="Progression-free survival,%",risk.table.title="",
legend.labs=c("SD/PD","CR/PR"),
legend.title="Subtype",risk.table.height = 0.3,surv.plot.height = 0.7
)
ggsurv$table <- ggsurv$table + theme_cleantable()
ggsurv

dev.off()




##############################################################################
## The sankey diagram between irAEs and response and pfs 

library(ggalluvial)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

clin_irAEs<-subset(clin,irAEs!="NA",drop=T)

nr<-nrow(clin_irAEs)
for(i in 1:nr){
  if(is.na(clin_irAEs$response_code[i])==T){clin_irAEs$clinical_resp_num[i]="Unknown"}
   else{
    if (clin_irAEs$response_code[i]==0){clin_irAEs$clinical_resp_num[i]="SD/PD"}
    if (clin_irAEs$response_code[i]==1){clin_irAEs$clinical_resp_num[i]="CR/PR"}
  }
}

nr<-nrow(clin_irAEs)
for(i in 1:nr){
  if(is.na(clin_irAEs$irAEs_grade[i])==T){clin_irAEs$Most_Severe_irAEs[i]="No irAEs"}
   else{
    #if (clin_irAEs$irAEs_grade[i]=="No_AE_reported"){clin_irAEs$Most_Severe_irAEs[i]="No irAEs"}
    if (clin_irAEs$irAEs_grade[i]=="G1-G2"){clin_irAEs$Most_Severe_irAEs[i]="Grade 1-2"}
    if (clin_irAEs$irAEs_grade[i]=="G3-G4"){clin_irAEs$Most_Severe_irAEs[i]="Grade 3-4"}
  }
}


nr<-nrow(clin_irAEs)
for(i in 1:nr){
  if(is.na(clin_irAEs$pfs_12_months[i])==T){clin_irAEs$pfs_12_months_nam[i]="Unknown"}
   else{
    if (clin_irAEs$pfs_12_months[i]==1){clin_irAEs$pfs_12_months_nam[i]="PFS >=12 months"}
    if (clin_irAEs$pfs_12_months[i]==0){clin_irAEs$pfs_12_months_nam[i]="PFS < 12 months"}
  }
}


head(clin_irAEs)
df<-clin_irAEs[,c("Most_Severe_irAEs","clinical_resp_num","pfs_12_months_nam")]
head(df)

#define color 
mycol<-c("#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#BC3C29FF","#EE4C97FF","#0072B5FF")

## Discrete Distribution Plot:
plot(table(rpois(100, 5)), type = "h", col = "#91612D", lwd = 10,
     main = "rpois(100, lambda = 5)")

#格式转换
UCB_lodes <- to_lodes_form(df[,1:ncol(df)],
                           axes = 1:ncol(df),
                           id = "Cohort")
dim(UCB_lodes)
head(UCB_lodes)
tail(UCB_lodes)

tiff(file = "pics/F1_ggalluvial_cluster.tiff", width =3500, height = 3000, res =600) 

ggplot(UCB_lodes,
       aes(x = x, stratum = stratum, alluvium = Cohort,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/6) + 
  geom_stratum(alpha = 1,width = 1/1.6,color="white") +
  geom_text(stat = "stratum", size =4,color="white") +   
  scale_fill_manual(values = mycol) +

  xlab("") + ylab("") +
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + 
  ggtitle("")+
  guides(fill = FALSE) 

dev.off()


##rev(mycolor6[-6])



################################################################################################
### association between clinical information and prognosis 

France_basic<-read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")
USA_UK_basic<-read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv")
clin<-rbind(France_basic,USA_UK_basic)

table(clin$cancer_type)
clin_resp<-subset(clin,response_code!="NA",drop=T)

lg_res_model <- glm(response_code~age,family=binomial(logit), data=clin_resp)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~gender,family=binomial(logit), data=clin_resp)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~BMI,family=binomial(logit), data=clin_resp)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~LDH,family=binomial(logit), data=clin_resp)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~cancer_type,family=binomial(logit), data=clin_resp)
lg_res <- summary(lg_res_model)


### different cancer type 

mel<-subset(clin_resp,cancer_type=="melanoma",drop=T)
ccr<-subset(clin_resp,cancer_type=="RCC",drop=T)
NSCLC<-subset(clin_resp,cancer_type=="NSCLC",drop=T)

lg_res_model <- glm(response_code~age,family=binomial(logit), data=mel)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~gender,family=binomial(logit), data=mel)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~BMI,family=binomial(logit), data=mel)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~LDH,family=binomial(logit), data=mel)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~age,family=binomial(logit), data=ccr)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~gender,family=binomial(logit), data=ccr)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~age,family=binomial(logit), data=NSCLC)
lg_res <- summary(lg_res_model)

lg_res_model <- glm(response_code~gender,family=binomial(logit), data=NSCLC)
lg_res <- summary(lg_res_model)



################################################################## os 
mel<-subset(clin,cancer_type=="melanoma",drop=T)
ccr<-subset(clin,cancer_type=="RCC",drop=T)
NSCLC<-subset(clin,cancer_type=="NSCLC",drop=T)

cox_res_model <- coxph(Surv(os, os_event)~age,clin,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(os, os_event)~gender,clin,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(os, os_event)~BMI,clin,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(os, os_event)~cancer_type,clin,na.action=na.exclude)
cox_res <- summary(cox_res_model)


cox_res_model <- coxph(Surv(os, os_event)~age,mel,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(os, os_event)~gender,mel,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(os, os_event)~BMI,mel,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(os, os_event)~age,NSCLC,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(os, os_event)~gender,NSCLC,na.action=na.exclude)
cox_res <- summary(cox_res_model)


################################################################## pfs
mel<-subset(clin,cancer_type=="melanoma",drop=T)
ccr<-subset(clin,cancer_type=="RCC",drop=T)
NSCLC<-subset(clin,cancer_type=="NSCLC",drop=T)

cox_res_model <- coxph(Surv(pfs, pfs_event)~age,clin,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(pfs, pfs_event)~gender,clin,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(pfs, pfs_event)~BMI,clin,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(pfs, pfs_event)~cancer_type,clin,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(pfs, pfs_event)~age,mel,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(pfs, pfs_event)~gender,mel,na.action=na.exclude)
cox_res <- summary(cox_res_model)

cox_res_model <- coxph(Surv(pfs, pfs_event)~BMI,mel,na.action=na.exclude)
cox_res <- summary(cox_res_model)


