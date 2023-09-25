### SV data processing
### 2023-9-19
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")


PRJNA770295_dsgv<-read.csv("00.rawData/SV_datasets/PRJNA770295_dsgv.csv",header=T,row.names=1)
t_PRJNA770295_dsgv<-data.frame(t(PRJNA770295_dsgv))
t_PRJNA770295_dsgv$id<-row.names(t_PRJNA770295_dsgv)

PRJNA397906_dsgv<-read.csv("00.rawData/SV_datasets/PRJNA397906_dsgv.csv",header=T,row.names=1)
t_PRJNA397906_dsgv<-data.frame(t(PRJNA397906_dsgv))
t_PRJNA397906_dsgv$id<-row.names(t_PRJNA397906_dsgv)

PRJEB43119_dsgv<-read.csv("00.rawData/SV_datasets/PRJEB43119_dsgv.csv",header=T,row.names=1)
t_PRJEB43119_dsgv<-data.frame(t(PRJEB43119_dsgv))
t_PRJEB43119_dsgv$id<-row.names(t_PRJEB43119_dsgv)

PRJNA541981_dsgv<-read.csv("00.rawData/SV_datasets/PRJNA541981_dsgv.csv",header=T,row.names=1)
t_PRJNA541981_dsgv<-data.frame(t(PRJNA541981_dsgv))
t_PRJNA541981_dsgv$id<-row.names(t_PRJNA541981_dsgv)

PRJNA762360_dsgv<-read.csv("00.rawData/SV_datasets/PRJNA762360_dsgv.csv",header=T,row.names=1)
t_PRJNA762360_dsgv<-data.frame(t(PRJNA762360_dsgv))
t_PRJNA762360_dsgv$id<-row.names(t_PRJNA762360_dsgv)

USA_UK_dsgv<-merge(t_PRJNA770295_dsgv,t_PRJNA397906_dsgv,by.x="id",by.y="id",all.x=T,all.y=T)
USA_UK_dsgv<-merge(USA_UK_dsgv,t_PRJNA541981_dsgv,by.x="id",by.y="id",all.x=T,all.y=T)
USA_UK_dsgv<-merge(USA_UK_dsgv,t_PRJNA762360_dsgv,by.x="id",by.y="id",all.x=T,all.y=T)
USA_UK_dsgv<-merge(USA_UK_dsgv,t_PRJEB43119_dsgv,by.x="id",by.y="id",all.x=T,all.y=T)

row.names(USA_UK_dsgv)<-USA_UK_dsgv$id
t_USA_UK_dsgv<-t(USA_UK_dsgv[,-1])
write.csv(t_USA_UK_dsgv,"00.rawData/SV/USA_UK_dsgv.csv")



PRJNA770295_vsgv<-read.csv("00.rawData/SV_datasets/PRJNA770295_vsgv.csv",header=T,row.names=1)
t_PRJNA770295_vsgv<-data.frame(t(PRJNA770295_vsgv))
t_PRJNA770295_vsgv$id<-row.names(t_PRJNA770295_vsgv)

PRJNA397906_vsgv<-read.csv("00.rawData/SV_datasets/PRJNA397906_vsgv.csv",header=T,row.names=1)
t_PRJNA397906_vsgv<-data.frame(t(PRJNA397906_vsgv))
t_PRJNA397906_vsgv$id<-row.names(t_PRJNA397906_vsgv)

PRJEB43119_vsgv<-read.csv("00.rawData/SV_datasets/PRJEB43119_vsgv.csv",header=T,row.names=1)
t_PRJEB43119_vsgv<-data.frame(t(PRJEB43119_vsgv))
t_PRJEB43119_vsgv$id<-row.names(t_PRJEB43119_vsgv)

PRJNA541981_vsgv<-read.csv("00.rawData/SV_datasets/PRJNA541981_vsgv.csv",header=T,row.names=1)
t_PRJNA541981_vsgv<-data.frame(t(PRJNA541981_vsgv))
t_PRJNA541981_vsgv$id<-row.names(t_PRJNA541981_vsgv)

PRJNA762360_vsgv<-read.csv("00.rawData/SV_datasets/PRJNA762360_vsgv.csv",header=T,row.names=1)
t_PRJNA762360_vsgv<-data.frame(t(PRJNA762360_vsgv))
t_PRJNA762360_vsgv$id<-row.names(t_PRJNA762360_vsgv)

USA_UK_vsgv<-merge(t_PRJNA770295_vsgv,t_PRJNA397906_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)
USA_UK_vsgv<-merge(USA_UK_vsgv,t_PRJNA541981_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)
USA_UK_vsgv<-merge(USA_UK_vsgv,t_PRJNA762360_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)
USA_UK_vsgv<-merge(USA_UK_vsgv,t_PRJEB43119_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)

row.names(USA_UK_vsgv)<-USA_UK_vsgv$id
t_USA_UK_vsgv<-t(USA_UK_vsgv[,-1])
write.csv(t_USA_UK_vsgv,"00.rawData/SV/USA_UK_vsgv.csv")




#####################################################################################
## NSCLC 

PRJEB22863_dsgv<-read.csv("00.rawData/SV_datasets/PRJEB22863_dsgv.csv",header=T,row.names=1)
t_PRJEB22863_dsgv<-data.frame(t(PRJEB22863_dsgv))
t_PRJEB22863_dsgv$id<-row.names(t_PRJEB22863_dsgv)

PRJNA751792_dsgv<-read.csv("00.rawData/SV_datasets/PRJNA751792_dsgv.csv",header=T,row.names=1)
t_PRJNA751792_dsgv<-data.frame(t(PRJNA751792_dsgv))
t_PRJNA751792_dsgv$id<-row.names(t_PRJNA751792_dsgv)

France_dsgv<-merge(t_PRJEB22863_dsgv,t_PRJNA751792_dsgv,by.x="id",by.y="id",all.x=T,all.y=T)
row.names(France_dsgv)<-France_dsgv$id
t_France_dsgv<-t(France_dsgv[,-1])
write.csv(t_France_dsgv,"00.rawData/SV/France_dsgv.csv")


PRJEB22863_vsgv<-read.csv("00.rawData/SV_datasets/PRJEB22863_vsgv.csv",header=T,row.names=1)
t_PRJEB22863_vsgv<-data.frame(t(PRJEB22863_vsgv))
t_PRJEB22863_vsgv$id<-row.names(t_PRJEB22863_vsgv)

PRJNA751792_vsgv<-read.csv("00.rawData/SV_datasets/PRJNA751792_vsgv.csv",header=T,row.names=1)
t_PRJNA751792_vsgv<-data.frame(t(PRJNA751792_vsgv))
t_PRJNA751792_vsgv$id<-row.names(t_PRJNA751792_vsgv)

France_vsgv<-merge(t_PRJEB22863_vsgv,t_PRJNA751792_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)
row.names(France_vsgv)<-France_vsgv$id
t_France_vsgv<-t(France_vsgv[,-1])
write.csv(t_France_vsgv,"00.rawData/SV/France_vsgv.csv")





