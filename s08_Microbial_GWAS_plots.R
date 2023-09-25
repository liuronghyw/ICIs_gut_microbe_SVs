### Microbial GWA 
### 2023-9-20
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
#options(scipen = 200)

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)
library(ggthemes)

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


###################################################################################################
### 3 Results
### 3.1 Clean results
#### 3.1.1 vSV associations

## vsv
## merge result tables
## mel (response, os ,pfs)
## nsclc(response,os)
## rcc (response)

melanoma_vsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_resp_meta_pvalue_adjust.csv",header=T)
melanoma_vsv_resp_meta_pvalue.sig.edge<-subset(melanoma_vsv_resp_meta_pvalue.edge,pvalue_meta <= 0.01 
& het_p >= 0.05 & p_count>=2,drop=T)
melanoma_vsv_resp_meta_pvalue.sig.anno.edge<-left_join(melanoma_vsv_resp_meta_pvalue.sig.edge, USA_UK_vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_resp_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/melanoma_vsv_response_adjAbun.sig.anno.csv",row.names=F)

melanoma_vsv_os_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_os_pvalue_adjust.csv",header=T)
melanoma_vsv_os_meta_pvalue.sig.edge<-subset(melanoma_vsv_os_meta_pvalue.edge,fdr.p<= 0.1,drop=T)
melanoma_vsv_os_meta_pvalue.sig.anno.edge<-left_join(melanoma_vsv_os_meta_pvalue.sig.edge, USA_UK_vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_os_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/melanoma_vsv_os_adjAbun.sig.anno.csv",row.names=F)

melanoma_vsv_pfs_12_months_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_pfs_12_months_meta_pvalue_adjust.csv",header=T)
melanoma_vsv_pfs_12_months_meta_pvalue.sig.edge<-subset(melanoma_vsv_pfs_12_months_meta_pvalue.edge,pvalue_meta <= 0.01 
& het_p > 0.05 & p_count>=2 ,drop=T)
melanoma_vsv_pfs_12_months_meta_pvalue.sig.anno.edge<-left_join(melanoma_vsv_pfs_12_months_meta_pvalue.sig.edge, USA_UK_vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_pfs_12_months_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/melanoma_vsv_pfs_12_months_adjAbun.sig.anno.csv",row.names=F)

melanoma_vsv_irAEs_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_irAEs_pvalue_adjust.csv",header=T)
melanoma_vsv_irAEs_meta_pvalue.sig.edge<-subset(melanoma_vsv_irAEs_meta_pvalue.edge,fdr.p<= 0.1,drop=T)
melanoma_vsv_irAEs_meta_pvalue.sig.anno.edge<-left_join(melanoma_vsv_irAEs_meta_pvalue.sig.edge, USA_UK_vsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_vsv_irAEs_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/melanoma_vsv_irAEs_adjAbun.sig.anno.csv",row.names=F)

NSCLC_vsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/nsclc_vsv_resp_meta_pvalue_adjust.csv",header=T)
NSCLC_vsv_resp_meta_pvalue.sig.edge<-subset(NSCLC_vsv_resp_meta_pvalue.edge,pvalue_meta <= 0.01 
& het_p > 0.05 & p_count>=2,drop=T)
NSCLC_vsv_resp_meta_pvalue.sig.anno.edge<-left_join(NSCLC_vsv_resp_meta_pvalue.sig.edge, France_vsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_vsv_resp_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/NSCLC_vsv_response_adjAbun.sig.anno.csv",row.names=F)

NSCLC_vsv_os_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/nsclc_vsv_os_pvalue_adjust.csv",header=T)
NSCLC_vsv_os_meta_pvalue.sig.edge<-subset(NSCLC_vsv_os_meta_pvalue.edge,fdr.p<= 0.1,drop=T)
NSCLC_vsv_os_meta_pvalue.sig.anno.edge<-left_join(NSCLC_vsv_os_meta_pvalue.sig.edge, France_vsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_vsv_os_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/NSCLC_vsv_os_adjAbun.sig.anno.csv",row.names=F)

RCC_vsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/rcc_vsv_resp_pvalue_adjust.csv",header=T)
RCC_vsv_resp_meta_pvalue.sig.edge<-subset(RCC_vsv_resp_meta_pvalue.edge,fdr.p<= 0.11,drop=T)
RCC_vsv_resp_meta_pvalue.sig.anno.edge<-left_join(RCC_vsv_resp_meta_pvalue.sig.edge, France_vsv_info, by = c("X" = "SV_Name"))
write.csv(RCC_vsv_resp_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/RCC_vsv_response_adjAbun.sig.anno.csv",row.names=F)


######################################################################################################
### dsv

melanoma_dsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_resp_meta_pvalue_adjust.csv",header=T)
melanoma_dsv_resp_meta_pvalue.sig.edge<-subset(melanoma_dsv_resp_meta_pvalue.edge,pvalue_meta <= 0.01 & het_p > 0.05 
& p_count>=2,drop=T)
melanoma_dsv_resp_meta_pvalue.sig.anno.edge<-left_join(melanoma_dsv_resp_meta_pvalue.sig.edge, USA_UK_dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_resp_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/melanoma_dsv_response_adjAbun.sig.anno.csv",row.names=F)

melanoma_dsv_os_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_os_pvalue_adjust.csv",header=T)
melanoma_dsv_os_meta_pvalue.sig.edge<-subset(melanoma_dsv_os_meta_pvalue.edge,fdr.p<= 0.1,drop=T)
melanoma_dsv_os_meta_pvalue.sig.anno.edge<-left_join(melanoma_dsv_os_meta_pvalue.sig.edge, USA_UK_dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_os_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/melanoma_dsv_os_adjAbun.sig.anno.csv",row.names=F)

melanoma_dsv_pfs_12_months_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_pfs_12_months_meta_pvalue_adjust.csv",header=T)
melanoma_dsv_pfs_12_months_meta_pvalue.sig.edge<-subset(melanoma_dsv_pfs_12_months_meta_pvalue.edge,pvalue_meta <= 0.01 
& het_p > 0.05 & p_count>=2,drop=T)
melanoma_dsv_pfs_12_months_meta_pvalue.sig.anno.edge<-left_join(melanoma_dsv_pfs_12_months_meta_pvalue.sig.edge, USA_UK_dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_pfs_12_months_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/melanoma_dsv_pfs_12_months_adjAbun.sig.anno.csv",row.names=F)

melanoma_dsv_irAEs_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_irAEs_pvalue_adjust.csv",header=T)
melanoma_dsv_irAEs_meta_pvalue.sig.edge<-subset(melanoma_dsv_irAEs_meta_pvalue.edge,fdr.p<= 0.1,drop=T)
melanoma_dsv_irAEs_meta_pvalue.sig.anno.edge<-left_join(melanoma_dsv_irAEs_meta_pvalue.sig.edge, USA_UK_dsv_info, by = c("X" = "SV_Name"))
write.csv(melanoma_dsv_irAEs_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/melanoma_dsv_irAEs_adjAbun.sig.anno.csv",row.names=F)

NSCLC_dsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/nsclc_dsv_resp_meta_pvalue_adjust.csv",header=T)
NSCLC_dsv_resp_meta_pvalue.sig.edge<-subset(NSCLC_dsv_resp_meta_pvalue.edge,pvalue_meta <= 0.01
& het_p > 0.05 & p_count>=2,drop=T)
NSCLC_dsv_resp_meta_pvalue.sig.anno.edge<-left_join(NSCLC_dsv_resp_meta_pvalue.sig.edge, France_dsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_dsv_resp_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/NSCLC_dsv_response_adjAbun.sig.anno.csv",row.names=F)

NSCLC_dsv_os_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/nsclc_dsv_os_pvalue_adjust.csv",header=T)
NSCLC_dsv_os_meta_pvalue.sig.edge<-subset(NSCLC_dsv_os_meta_pvalue.edge,fdr.p<= 0.1,drop=T)
NSCLC_dsv_os_meta_pvalue.sig.anno.edge<-left_join(NSCLC_dsv_os_meta_pvalue.sig.edge, France_dsv_info, by = c("X" = "SV_Name"))
write.csv(NSCLC_dsv_os_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/NSCLC_dsv_os_adjAbun.sig.anno.csv",row.names=F)

RCC_dsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/rcc_dsv_resp_pvalue_adjust.csv",header=T)
RCC_dsv_resp_meta_pvalue.sig.edge<-subset(RCC_dsv_resp_meta_pvalue.edge,fdr.p<= 0.1,drop=T)
RCC_dsv_resp_meta_pvalue.sig.anno.edge<-left_join(RCC_dsv_resp_meta_pvalue.sig.edge, France_dsv_info, by = c("X" = "SV_Name"))
write.csv(RCC_dsv_resp_meta_pvalue.sig.anno.edge,"07.Microbial_GWAS/RCC_dsv_response_adjAbun.sig.anno.csv",row.names=F)


########################################################
#####  3.3 Visualization  （将不同的表型的结果汇总起来）
####   3.3.1 Combine vSV and dSV associations

melanoma_vsv_resp_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/melanoma_vsv_response_adjAbun.sig.anno.csv",header=T)
melanoma_vsv_os_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/melanoma_vsv_os_adjAbun.sig.anno.csv",header=T)
melanoma_vsv_pfs_12_months_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/melanoma_vsv_pfs_12_months_adjAbun.sig.anno.csv",header=T)
melanoma_vsv_irAEs_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/melanoma_vsv_irAEs_adjAbun.sig.anno.csv",header=T)

melanoma_dsv_resp_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/melanoma_dsv_response_adjAbun.sig.anno.csv",header=T)
melanoma_dsv_os_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/melanoma_dsv_os_adjAbun.sig.anno.csv",header=T)
melanoma_dsv_pfs_12_months_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/melanoma_dsv_pfs_12_months_adjAbun.sig.anno.csv",header=T)
melanoma_dsv_irAEs_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/melanoma_dsv_irAEs_adjAbun.sig.anno.csv",header=T)

melanoma_vsv_resp_meta_pvalue.sig.anno.edge$Pheno<-"Melanoma\nResponse"
#melanoma_vsv_os_meta_pvalue.sig.anno.edge$Pheno<-"Melanoma\nOS"
melanoma_vsv_pfs_12_months_meta_pvalue.sig.anno.edge$Pheno<-"Melanoma\nPFS >=12 months"
#melanoma_vsv_irAEs_meta_pvalue.sig.anno.edge$Pheno<-"Melanoma\nirAEs"

melanoma_dsv_resp_meta_pvalue.sig.anno.edge$Pheno<-"Melanoma\nResponse"
#melanoma_dsv_os_meta_pvalue.sig.anno.edge$Pheno<-"Melanoma\nOS"
melanoma_dsv_pfs_12_months_meta_pvalue.sig.anno.edge$Pheno<-"Melanoma\nPFS >=12 months"
#melanoma_dsv_irAEs_meta_pvalue.sig.anno.edge$Pheno<-"Melanoma\nirAEs"

melanoma_dsv_resp_adjAbun_sub<-melanoma_dsv_resp_meta_pvalue.sig.anno.edge[,c("Pheno","X","pvalue_meta","het_p")]
melanoma_dsv_pfs_12_months_adjAbun_sub<-melanoma_dsv_pfs_12_months_meta_pvalue.sig.anno.edge[,c("Pheno","X","pvalue_meta","het_p")]

melanoma_vsv_resp_adjAbun_sub<-melanoma_vsv_resp_meta_pvalue.sig.anno.edge[,c("Pheno","X","pvalue_meta","het_p")]
melanoma_vsv_pfs_12_months_adjAbun_sub<-melanoma_vsv_pfs_12_months_meta_pvalue.sig.anno.edge[,c("Pheno","X","pvalue_meta","het_p")]

#melanoma_vsv_os_adjAbun_sub<-melanoma_vsv_os_meta_pvalue.sig.anno.edge[,c("Pheno","X","fdr.p")]
#melanoma_vsv_os_adjAbun_sub$pvalue_meta<-melanoma_vsv_os_adjAbun_sub$fdr.p
#melanoma_vsv_os_adjAbun_sub$het_p<-1
#melanoma_vsv_os_adjAbun_sub<-melanoma_vsv_os_adjAbun_sub[,c("Pheno","X","pvalue_meta","het_p")]

#melanoma_dsv_os_adjAbun_sub<-melanoma_dsv_os_meta_pvalue.sig.anno.edge[,c("Pheno","X","fdr.p")]
#melanoma_dsv_os_adjAbun_sub$pvalue_meta<-melanoma_dsv_os_adjAbun_sub$fdr.p
#melanoma_dsv_os_adjAbun_sub$het_p<-1
#melanoma_dsv_os_adjAbun_sub<-melanoma_dsv_os_adjAbun_sub[,c("Pheno","X","pvalue_meta","het_p")]

#melanoma_vsv_irAEs_adjAbun_sub<-melanoma_vsv_irAEs_meta_pvalue.sig.anno.edge[,c("Pheno","X","p")]
#melanoma_vsv_irAEs_adjAbun_sub$pvalue_meta<-melanoma_vsv_irAEs_adjAbun_sub$p
#melanoma_vsv_irAEs_adjAbun_sub$het_p<-1
#melanoma_vsv_irAEs_adjAbun_sub<-melanoma_vsv_irAEs_adjAbun_sub[,c("Pheno","X","pvalue_meta","het_p")]

#melanoma_dsv_irAEs_adjAbun_sub<-melanoma_dsv_irAEs_meta_pvalue.sig.anno.edge[,c("Pheno","X","p")]
#melanoma_dsv_irAEs_adjAbun_sub$pvalue_meta<-melanoma_dsv_irAEs_adjAbun_sub$p
#melanoma_dsv_irAEs_adjAbun_sub$het_p<-1
#melanoma_dsv_irAEs_adjAbun_sub<-melanoma_dsv_irAEs_adjAbun_sub[,c("Pheno","X","pvalue_meta","het_p")]

#sv_prog_meta_pvalue.sig.edge<-rbind(melanoma_vsv_os_adjAbun_sub,melanoma_dsv_os_adjAbun_sub,
#melanoma_vsv_resp_adjAbun_sub,melanoma_dsv_resp_adjAbun_sub,
#melanoma_vsv_pfs_12_months_adjAbun_sub,melanoma_dsv_pfs_12_months_adjAbun_sub)

sv_prog_meta_pvalue.sig.edge<-rbind(melanoma_dsv_resp_adjAbun_sub,melanoma_vsv_resp_adjAbun_sub,
melanoma_dsv_pfs_12_months_adjAbun_sub,melanoma_vsv_pfs_12_months_adjAbun_sub)


colnames(sv_prog_meta_pvalue.sig.edge)[2]<-"SV"
sv_prog_adjAbun.sig.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno %>%
  as.character(.) %>%
  .[!duplicated(.)]

sv_prog_adjAbun.sig.sv <- sv_prog_meta_pvalue.sig.edge$SV %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

## circos plot
sv_prog_count<-NULL

for (phe in sv_prog_adjAbun.sig.pheno) {
  #pheno <-"C4"
  sv_prog_df<-sv_prog_meta_pvalue.sig.edge[sv_prog_meta_pvalue.sig.edge$Pheno==phe,]$SV %>%
    str_replace_all(.,"\\:\\d+_\\d+.*", "") %>%
    table %>%
    as.data.frame
  colnames(sv_prog_df)<-c("Species", "Count")
  sv_prog_df<-data.frame(Pheno = rep(phe,nrow(sv_prog_df)), sv_prog_df)
  sv_prog_count<-rbind(sv_prog_count,sv_prog_df)
}

sv_prog_count$Species<-USA_UK_info$Short_name[match(sv_prog_count$Species, USA_UK_info$organism)]
sv_prog_count <- sv_prog_count[order(sv_prog_count$Count),]

sv_prog_count_species_order<-sv_prog_count %>% group_by(Species) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]
sv_prog_count_pheno_order<-sv_prog_count %>% group_by(Pheno) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]

sv_prog_count<-sv_prog_count[order(match(sv_prog_count$Pheno, sv_prog_count_pheno_order$Pheno),decreasing = F),]

sv_prog_count_species_order_str<-as.character(sv_prog_count_species_order$Species)
sv_prog_count_pheno_order_str<-as.character(sv_prog_count_pheno_order$Pheno)


#pdf("pics/F6_SV_pheno_adjAbun.circos.pdf", width = 26, height = 26)
tiff(file = "pics/F5_SV_pheno_Mel_adjAbun.circos.tiff", width =7000, height =6300, res =300) 

circos.clear()
circos.par(start.degree = 180, "clock.wise" = T) # ,xaxis.clock.wise = F

chordDiagram(sv_prog_count,annotationTrack = "grid",
             #grid.col =  c(wes_palette("Darjeeling1", length(sv_prog_count_pheno_order_str), type = "continuous"),
             #              rep('grey',length(sv_prog_count_species_order_str))),
             grid.col =  c(c("#0072B5FF","#EE4C97FF"),
                           rep('grey',length(sv_prog_count_species_order_str))),
             order = c(rev(sv_prog_count_pheno_order_str),
                       sv_prog_count_species_order_str),
             big.gap = 6,
             preAllocateTracks = list(track.margin = c(0, uh(45, "mm")), 
             #                         #track.height = max(strwidth(unlist(dimnames(sv_prog_count)))))
                                       track.height=0.1)
             )

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,cex = 2.5,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

dev.off()


#mycolor6
hist(mtcars$mpg,
     breaks = 15,
     col = "#6F99ADFF",
)

## "#7876B1FF","#EE4C97FF"
#"#6F99ADFF","#E18727FF","#BC3C29FF","#7876B1FF","#0072B5FF",#20854EFF,"#EE4C97FF"

write.csv(sv_prog_count,"07.Microbial_GWAS/mel_sv_prog_count.csv",row.names=F)
unique(sv_prog_count$Species)
sum(sv_prog_count$Count)

unique(sv_prog_count$Species)
unique(sv_prog_meta_pvalue.sig.edge$SV)


############# 3.4 Examples
### 3.4.1 Heatmap of certain species

melanoma_vsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_resp_meta_pvalue.csv",header=T)
melanoma_vsv_resp_meta_pvalue.anno.edge<-left_join(melanoma_vsv_resp_meta_pvalue.edge, USA_UK_vsv_info, by = c("X" = "SV_Name"))

melanoma_vsv_os_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_os_pvalue_adjust.csv",header=T)
melanoma_vsv_os_meta_pvalue.anno.edge<-left_join(melanoma_vsv_os_meta_pvalue.edge, USA_UK_vsv_info, by = c("X" = "SV_Name"))

melanoma_vsv_pfs_12_months_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_pfs_12_months_meta_pvalue.csv",header=T)
melanoma_vsv_pfs_12_months_meta_pvalue.anno.edge<-left_join(melanoma_vsv_pfs_12_months_meta_pvalue.edge, USA_UK_vsv_info, by = c("X" = "SV_Name"))

melanoma_vsv_resp_adjAbun_sig_spe<-subset(melanoma_vsv_resp_meta_pvalue.edge,as.numeric(pvalue_meta) <= 0.01 & het_p > 0.05 & p_count >=2,drop=T)$X
melanoma_vsv_os_adjAbun_sig_spe<-subset(melanoma_vsv_os_meta_pvalue.edge,as.numeric(fdr.p) <= 0.1,drop=T)$X
melanoma_vsv_pfs_12_months_adjAbun_sig_spe<-subset(melanoma_vsv_pfs_12_months_meta_pvalue.edge,as.numeric(pvalue_meta) < 0.01 & het_p > 0.05 & p_count >=2,drop=T)$X

vsv_select<-data.frame(c(melanoma_vsv_os_adjAbun_sig_spe,melanoma_vsv_resp_adjAbun_sig_spe,melanoma_vsv_pfs_12_months_adjAbun_sig_spe))
colnames(vsv_select)<-"X"

melanoma_dsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_resp_meta_pvalue.csv",header=T)
melanoma_dsv_resp_meta_pvalue.anno.edge<-left_join(melanoma_dsv_resp_meta_pvalue.edge, USA_UK_dsv_info, by = c("X" = "SV_Name"))

melanoma_dsv_os_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_os_pvalue_adjust.csv",header=T)
melanoma_dsv_os_meta_pvalue.anno.edge<-left_join(melanoma_dsv_os_meta_pvalue.edge, USA_UK_dsv_info, by = c("X" = "SV_Name"))

melanoma_dsv_pfs_12_months_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/mel_dsv_pfs_12_months_meta_pvalue.csv",header=T)
melanoma_dsv_pfs_12_months_meta_pvalue.anno.edge<-left_join(melanoma_dsv_pfs_12_months_meta_pvalue.edge, USA_UK_dsv_info, by = c("X" = "SV_Name"))

melanoma_dsv_resp_adjAbun_sig_spe<-subset(melanoma_dsv_resp_meta_pvalue.edge,as.numeric(pvalue_meta) < 0.01 & het_p > 0.05 & p_count >=2,drop=T)$X
melanoma_dsv_os_adjAbun_sig_spe<-subset(melanoma_dsv_os_meta_pvalue.edge,as.numeric(fdr.p) <= 0.1,drop=T)$X
melanoma_dsv_pfs_12_months_adjAbun_sig_spe<-subset(melanoma_dsv_pfs_12_months_meta_pvalue.edge,as.numeric(pvalue_meta) < 0.01 & het_p > 0.05 & p_count >=2,drop=T)$X

dsv_select<-data.frame(c(melanoma_dsv_os_adjAbun_sig_spe,melanoma_dsv_resp_adjAbun_sig_spe,melanoma_dsv_pfs_12_months_adjAbun_sig_spe))
colnames(dsv_select)<-"X"

melanoma_vsv_resp_adjAbun_select<-merge(vsv_select,melanoma_vsv_resp_meta_pvalue.anno.edge,by.x="X",by.y="X")
melanoma_vsv_pfs_12_months_adjAbun_select<-merge(vsv_select,melanoma_vsv_pfs_12_months_meta_pvalue.anno.edge,by.x="X",by.y="X")
melanoma_vsv_os_adjAbun_select<-merge(vsv_select,melanoma_vsv_os_meta_pvalue.anno.edge,by.x="X",by.y="X")

melanoma_dsv_resp_adjAbun_select<-merge(dsv_select,melanoma_dsv_resp_meta_pvalue.anno.edge,by.x="X",by.y="X")
melanoma_dsv_pfs_12_months_adjAbun_select<-merge(dsv_select,melanoma_dsv_pfs_12_months_meta_pvalue.anno.edge,by.x="X",by.y="X")
melanoma_dsv_os_adjAbun_select<-merge(dsv_select,melanoma_dsv_os_meta_pvalue.anno.edge,by.x="X",by.y="X")


# Melanoma.
melanoma_vsv_resp_adjAbun_select$Pheno<-"Response"
melanoma_vsv_os_adjAbun_select$Pheno<-"OS"
melanoma_vsv_pfs_12_months_adjAbun_select$Pheno<-"PFS"

melanoma_dsv_resp_adjAbun_select$Pheno<-"Response"
melanoma_dsv_os_adjAbun_select$Pheno<-"OS"
melanoma_dsv_pfs_12_months_adjAbun_select$Pheno<-"PFS"

melanoma_vsv_resp_adjAbun_sub<-melanoma_vsv_resp_adjAbun_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]
melanoma_vsv_pfs_12_months_adjAbun_sub<-melanoma_vsv_pfs_12_months_adjAbun_select[,c("Pheno","X","pvalue_meta","se_meta","effect_meta","cilb","ciub")]

melanoma_vsv_os_adjAbun_sub<-melanoma_vsv_os_adjAbun_select[,c("Pheno","X","fdr.p","HR","SE","HR_left","HR_right")]
melanoma_vsv_os_adjAbun_sub$pvalue_meta<-melanoma_vsv_os_adjAbun_sub$fdr.p
melanoma_vsv_os_adjAbun_sub$het_p<-1
melanoma_vsv_os_adjAbun_sub$effect_meta<-melanoma_vsv_os_adjAbun_sub$HR
melanoma_vsv_os_adjAbun_sub$cilb<-melanoma_vsv_os_adjAbun_sub$HR_left
melanoma_vsv_os_adjAbun_sub$ciub<-melanoma_vsv_os_adjAbun_sub$HR_right
melanoma_vsv_os_adjAbun_sub$se_meta<-melanoma_vsv_os_adjAbun_sub$SE
melanoma_vsv_os_adjAbun_sub<-melanoma_vsv_os_adjAbun_sub[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]

melanoma_dsv_resp_adjAbun_sub<-melanoma_dsv_resp_adjAbun_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]
melanoma_dsv_pfs_12_months_adjAbun_sub<-melanoma_dsv_pfs_12_months_adjAbun_select[,c("Pheno","X","pvalue_meta","se_meta","effect_meta","cilb","ciub")]

melanoma_dsv_os_adjAbun_sub<-melanoma_dsv_os_adjAbun_select[,c("Pheno","X","fdr.p","HR","SE","HR_left","HR_right")]
melanoma_dsv_os_adjAbun_sub$pvalue_meta<-melanoma_dsv_os_adjAbun_sub$fdr.p
melanoma_dsv_os_adjAbun_sub$het_p<-1
melanoma_dsv_os_adjAbun_sub$effect_meta<-melanoma_dsv_os_adjAbun_sub$HR
melanoma_dsv_os_adjAbun_sub$cilb<-melanoma_dsv_os_adjAbun_sub$HR_left
melanoma_dsv_os_adjAbun_sub$ciub<-melanoma_dsv_os_adjAbun_sub$HR_right
melanoma_dsv_os_adjAbun_sub$se_meta<-melanoma_dsv_os_adjAbun_sub$SE
melanoma_dsv_os_adjAbun_sub<-melanoma_dsv_os_adjAbun_sub[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]

melanoma_vsv_resp_adjAbun_sub$SV<-c(paste("vsv:",melanoma_vsv_resp_adjAbun_sub$X,sep = ""))
melanoma_vsv_os_adjAbun_sub$SV<-c(paste("vsv:",melanoma_vsv_os_adjAbun_sub$X,sep = ""))
melanoma_vsv_pfs_12_months_adjAbun_sub$SV<-c(paste("vsv:",melanoma_vsv_pfs_12_months_adjAbun_sub$X,sep = ""))

melanoma_dsv_resp_adjAbun_sub$SV<-c(paste("dsv:",melanoma_dsv_resp_adjAbun_sub$X,sep = ""))
melanoma_dsv_os_adjAbun_sub$SV<-c(paste("dsv:",melanoma_dsv_os_adjAbun_sub$X,sep = ""))
melanoma_dsv_pfs_12_months_adjAbun_sub$SV<-c(paste("dsv:",melanoma_dsv_pfs_12_months_adjAbun_sub$X,sep = ""))

sv_prog_meta_pvalue.sig.edge<-rbind(melanoma_vsv_resp_adjAbun_sub,
melanoma_vsv_pfs_12_months_adjAbun_sub,melanoma_dsv_resp_adjAbun_sub,
melanoma_dsv_pfs_12_months_adjAbun_sub,melanoma_dsv_os_adjAbun_sub,melanoma_vsv_os_adjAbun_sub)

sv_prog_meta_pvalue.sig.edge$beta<-log(sv_prog_meta_pvalue.sig.edge$effect_meta)
sv_prog_meta_pvalue.sig.edge$beta_left<-log(sv_prog_meta_pvalue.sig.edge$cilb)
sv_prog_meta_pvalue.sig.edge$beta_right<-log(sv_prog_meta_pvalue.sig.edge$ciub)

sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge$Sign<-NA

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.1 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "*" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.05 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "☆" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.01 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "★" 

sv_prog_meta_pvalue.sig.edge<-unique(sv_prog_meta_pvalue.sig.edge)

sv_prog_adjAbun.r<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "beta")] %>% spread("Pheno", "beta") %>% data.frame(row.names = "SV")
sv_prog_adjAbun.sign<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "Sign")] %>% spread("Pheno", "Sign") %>% data.frame(row.names = "SV")

sv_prog_meta_pvalue.sig.edge<-subset(sv_prog_meta_pvalue.sig.edge,is.na(pvalue_meta)==F,drop=T)
# write.csv(sv_prog_meta_pvalue.sig.edge,"test.csv")


##########################################################################################################
###  A.muciniphila

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("muciniphila", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("muciniphila", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("muciniphila", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("muciniphila", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Akkermansia muciniphila ATCC BAA-835:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

#sv_adjAbun.r.plot<-sv_adjAbun.r.plot[,-3]

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-2, 0, length.out=ceiling(i/2) + 1), 
              seq(1/i, 2, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F5B_A.muciniphila_sv_adjAbun.heatmap.tiff", width =4200, height =2600, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("Response","PFS >= 12 months"),
          cexCol =3.5, srtCol = 45, cexRow = 3.5,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =5,
          key.xlab = "Beta coefficient",key=F,
          key.par=list(mar=c(4,12,4,4), cex.axis = 3.5, cex.lab = 3.5), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0)), lhei = c(0.5, 1.5),lwid=c(0.5, 1.5, 1 ),key.title = NA,
          margins=c(20,5) # ("margin.Y", "margin.X")
          )

dev.off()


##########################################################################################################
###  P.distasonis

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("distasonis", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("distasonis", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("distasonis", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("distasonis", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Parabacteroides distasonis ATCC 8503:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-2, 0, length.out=ceiling(i/2) + 1), 
              seq(2/i, 2, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

b<-b[,-3]

tiff(file = "pics/F5B_Pdistasonis_sv_adjAbun.heatmap.tiff", width =6000, height =1800, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("Response","PFS >= 12 months","OS"),
          cexCol =3.5, srtCol = 45, cexRow = 3.5,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =6,
          key.xlab = "Beta coefficient",key=F,
          key.par=list(mar=c(4,6,4,4), cex.axis = 3.5, cex.lab = 3.5), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 1.5, 1 ),key.title = NA,
          margins=c(20,20) # ("margin.Y", "margin.X")
          )

dev.off()



###############################################################################################
## NSCLC & RCC

NSCLC_vsv_resp_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/NSCLC_vsv_response_adjAbun.sig.anno.csv",header=T)
NSCLC_vsv_os_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/NSCLC_vsv_os_adjAbun.sig.anno.csv",header=T)

NSCLC_dsv_resp_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/NSCLC_dsv_response_adjAbun.sig.anno.csv",header=T)
NSCLC_dsv_os_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/NSCLC_dsv_os_adjAbun.sig.anno.csv",header=T)

RCC_vsv_resp_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/RCC_vsv_response_adjAbun.sig.anno.csv",header=T)
RCC_dsv_resp_meta_pvalue.sig.anno.edge<-read.csv("07.Microbial_GWAS/RCC_dsv_response_adjAbun.sig.anno.csv",header=T)

NSCLC_vsv_resp_meta_pvalue.sig.anno.edge$Pheno<-"NSCLC\nResponse"
NSCLC_vsv_os_meta_pvalue.sig.anno.edge$Pheno<-"NSCLC\nOS"

NSCLC_dsv_resp_meta_pvalue.sig.anno.edge$Pheno<-"NSCLC\nResponse"
NSCLC_dsv_os_meta_pvalue.sig.anno.edge$Pheno<-"NSCLC\nOS"

RCC_dsv_resp_meta_pvalue.sig.anno.edge$Pheno<-"RCC\nResponse"
RCC_vsv_resp_meta_pvalue.sig.anno.edge$Pheno<-"RCC\nResponse"

NSCLC_vsv_resp_adjAbun_sub<-NSCLC_vsv_resp_meta_pvalue.sig.anno.edge[,c("Pheno","X","pvalue_meta","het_p")]

NSCLC_vsv_os_adjAbun_sub<-NSCLC_vsv_os_meta_pvalue.sig.anno.edge[,c("Pheno","X","fdr.p")]
NSCLC_vsv_os_adjAbun_sub$pvalue_meta<-NSCLC_vsv_os_adjAbun_sub$fdr.p
NSCLC_vsv_os_adjAbun_sub$het_p<-1
NSCLC_vsv_os_adjAbun_sub<-NSCLC_vsv_os_adjAbun_sub[,c("Pheno","X","pvalue_meta","het_p")]

NSCLC_dsv_resp_adjAbun_sub<-NSCLC_dsv_resp_meta_pvalue.sig.anno.edge[,c("Pheno","X","pvalue_meta","het_p")]
NSCLC_dsv_resp_adjAbun_sub<-subset(NSCLC_dsv_resp_adjAbun_sub,is.na(pvalue_meta)==F,drop=T)

NSCLC_dsv_os_adjAbun_sub<-NSCLC_dsv_os_meta_pvalue.sig.anno.edge[,c("Pheno","X","fdr.p")]
NSCLC_dsv_os_adjAbun_sub$pvalue_meta<-NSCLC_dsv_os_adjAbun_sub$fdr.p
NSCLC_dsv_os_adjAbun_sub$het_p<-1
NSCLC_dsv_os_adjAbun_sub<-NSCLC_dsv_os_adjAbun_sub[,c("Pheno","X","pvalue_meta","het_p")]

RCC_dsv_resp_adjAbun_sub<-RCC_dsv_resp_meta_pvalue.sig.anno.edge[,c("Pheno","X","fdr.p")]
RCC_dsv_resp_adjAbun_sub$pvalue_meta<-RCC_dsv_resp_adjAbun_sub$fdr.p
RCC_dsv_resp_adjAbun_sub$het_p<-1
RCC_dsv_resp_adjAbun_sub<-RCC_dsv_resp_adjAbun_sub[,c("Pheno","X","pvalue_meta","het_p")]

RCC_vsv_resp_adjAbun_sub<-RCC_vsv_resp_meta_pvalue.sig.anno.edge[,c("Pheno","X","fdr.p")]
RCC_vsv_resp_adjAbun_sub$pvalue_meta<-RCC_vsv_resp_adjAbun_sub$fdr.p
RCC_vsv_resp_adjAbun_sub$het_p<-1
RCC_vsv_resp_adjAbun_sub<-RCC_vsv_resp_adjAbun_sub[,c("Pheno","X","pvalue_meta","het_p")]

sv_prog_meta_pvalue.sig.edge<-rbind(NSCLC_dsv_resp_adjAbun_sub,NSCLC_vsv_resp_adjAbun_sub,
NSCLC_vsv_os_adjAbun_sub,NSCLC_dsv_os_adjAbun_sub,RCC_dsv_resp_adjAbun_sub,RCC_vsv_resp_adjAbun_sub)

sv_prog_meta_pvalue.sig.edge<-subset(sv_prog_meta_pvalue.sig.edge,pvalue_meta<=0.1,drop=T)


colnames(sv_prog_meta_pvalue.sig.edge)[2]<-"SV"

sv_prog_adjAbun.sig.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno %>%
  as.character(.) %>%
  .[!duplicated(.)]

sv_prog_adjAbun.sig.sv <- sv_prog_meta_pvalue.sig.edge$SV %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

## circos plot
sv_prog_count<-NULL

for (phe in sv_prog_adjAbun.sig.pheno) {
  #pheno <-"C4"
  sv_prog_df<-sv_prog_meta_pvalue.sig.edge[sv_prog_meta_pvalue.sig.edge$Pheno==phe,]$SV %>%
    str_replace_all(.,"\\:\\d+_\\d+.*", "") %>%
    table %>%
    as.data.frame
  colnames(sv_prog_df)<-c("Species", "Count")
  sv_prog_df<-data.frame(Pheno = rep(phe,nrow(sv_prog_df)), sv_prog_df)
  sv_prog_count<-rbind(sv_prog_count,sv_prog_df)
}

sv_prog_count$Species<-France_info$Short_name[match(sv_prog_count$Species, France_info$organism)]
sv_prog_count <- sv_prog_count[order(sv_prog_count$Count),]

sv_prog_count_species_order<-sv_prog_count %>% group_by(Species) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]
sv_prog_count_pheno_order<-sv_prog_count %>% group_by(Pheno) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]

sv_prog_count<-sv_prog_count[order(match(sv_prog_count$Pheno, sv_prog_count_pheno_order$Pheno),decreasing = F),]

sv_prog_count_species_order_str<-as.character(sv_prog_count_species_order$Species)
sv_prog_count_pheno_order_str<-as.character(sv_prog_count_pheno_order$Pheno)


#pdf("pics/F6_SV_pheno_adjAbun.circos.pdf", width = 26, height = 26)
tiff(file = "pics/F6_SV_NSCLC_pheno_adjAbun.circos.tiff", width =7000, height =5800, res =300) 

circos.clear()
circos.par(start.degree = 180, "clock.wise" = T) # ,xaxis.clock.wise = F

chordDiagram(sv_prog_count,annotationTrack = "grid",
             #grid.col =  c(wes_palette("Darjeeling1", length(sv_prog_count_pheno_order_str), type = "continuous"),
             #              rep('grey',length(sv_prog_count_species_order_str))),
             grid.col =  c(c("#20854EFF","#0072B5FF","#7876B1FF"),
                           rep('grey',length(sv_prog_count_species_order_str))),
             order = c(rev(sv_prog_count_pheno_order_str),
                       sv_prog_count_species_order_str),
             big.gap = 6,
             preAllocateTracks = list(track.margin = c(0, uh(45, "mm")), 
             #                         #track.height = max(strwidth(unlist(dimnames(sv_prog_count)))))
                                       track.height=0.1)
             )

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,cex = 2.5,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

dev.off()

#mycolor6
#hist(mtcars$mpg,
#     breaks = 15,
#     col = "#EE4C97FF",
#)
#"#6F99ADFF","#E18727FF","#BC3C29FF","#7876B1FF","#0072B5FF",#20854EFF

write.csv(sv_prog_count,"07.Microbial_GWAS/NSCLC_sv_prog_count.csv",row.names=F)
unique(sv_prog_count$Species)
unique(sv_prog_meta_pvalue.sig.edge$SV)


### 3.4 Examples
### 3.4.1 Heatmap of certain species

NSCLC_vsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/NSCLC_vsv_resp_meta_pvalue.csv",header=T)
NSCLC_vsv_resp_meta_pvalue.anno.edge<-left_join(NSCLC_vsv_resp_meta_pvalue.edge, France_vsv_info, by = c("X" = "SV_Name"))

NSCLC_vsv_os_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/nsclc_vsv_os_pvalue_adjust.csv",header=T)
NSCLC_vsv_os_meta_pvalue.anno.edge<-left_join(NSCLC_vsv_os_meta_pvalue.edge, France_vsv_info, by = c("X" = "SV_Name"))

NSCLC_vsv_resp_adjAbun_sig_spe<-subset(NSCLC_vsv_resp_meta_pvalue.edge,as.numeric(pvalue_meta) <= 0.01 & het_p > 0.05 & p_count >=2,drop=T)$X
NSCLC_vsv_os_adjAbun_sig_spe<-subset(NSCLC_vsv_os_meta_pvalue.edge,as.numeric(fdr.p) <= 0.1,drop=T)$X

vsv_select<-data.frame(c(NSCLC_vsv_resp_adjAbun_sig_spe,NSCLC_vsv_os_adjAbun_sig_spe))
colnames(vsv_select)<-"X"

NSCLC_dsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/NSCLC_dsv_resp_meta_pvalue.csv",header=T)
NSCLC_dsv_resp_meta_pvalue.anno.edge<-left_join(NSCLC_dsv_resp_meta_pvalue.edge, France_dsv_info, by = c("X" = "SV_Name"))

NSCLC_dsv_os_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/nsclc_dsv_os_pvalue_adjust.csv",header=T)
NSCLC_dsv_os_meta_pvalue.anno.edge<-left_join(NSCLC_dsv_os_meta_pvalue.edge, France_dsv_info, by = c("X" = "SV_Name"))

NSCLC_dsv_resp_adjAbun_sig_spe<-subset(NSCLC_dsv_resp_meta_pvalue.edge,as.numeric(pvalue_meta) < 0.01 & het_p > 0.05 & p_count >= 2,drop=T)$X
NSCLC_dsv_os_adjAbun_sig_spe<-subset(NSCLC_dsv_os_meta_pvalue.edge,as.numeric(fdr.p) <= 0.1,drop=T)$X

dsv_select<-data.frame(c(NSCLC_dsv_resp_adjAbun_sig_spe,NSCLC_dsv_os_adjAbun_sig_spe))
colnames(dsv_select)<-"X"

NSCLC_vsv_resp_adjAbun_select<-merge(vsv_select,NSCLC_vsv_resp_meta_pvalue.anno.edge,by.x="X",by.y="X")
NSCLC_vsv_os_adjAbun_select<-merge(vsv_select,NSCLC_vsv_os_meta_pvalue.anno.edge,by.x="X",by.y="X")

NSCLC_dsv_resp_adjAbun_select<-merge(dsv_select,NSCLC_dsv_resp_meta_pvalue.anno.edge,by.x="X",by.y="X")
NSCLC_dsv_os_adjAbun_select<-merge(dsv_select,NSCLC_dsv_os_meta_pvalue.anno.edge,by.x="X",by.y="X")

NSCLC_vsv_resp_adjAbun_select$Pheno<-"Response"
NSCLC_vsv_os_adjAbun_select$Pheno<-"OS"

NSCLC_dsv_resp_adjAbun_select$Pheno<-"Response"
NSCLC_dsv_os_adjAbun_select$Pheno<-"OS"

NSCLC_vsv_resp_adjAbun_sub<-NSCLC_vsv_resp_adjAbun_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]

NSCLC_vsv_os_adjAbun_sub<-NSCLC_vsv_os_adjAbun_select[,c("Pheno","X","fdr.p","HR","SE","HR_left","HR_right")]
NSCLC_vsv_os_adjAbun_sub$pvalue_meta<-NSCLC_vsv_os_adjAbun_sub$fdr.p
NSCLC_vsv_os_adjAbun_sub$het_p<-1
NSCLC_vsv_os_adjAbun_sub$effect_meta<-NSCLC_vsv_os_adjAbun_sub$HR
NSCLC_vsv_os_adjAbun_sub$cilb<-NSCLC_vsv_os_adjAbun_sub$HR_left
NSCLC_vsv_os_adjAbun_sub$ciub<-NSCLC_vsv_os_adjAbun_sub$HR_right
NSCLC_vsv_os_adjAbun_sub$se_meta<-NSCLC_vsv_os_adjAbun_sub$SE

NSCLC_vsv_os_adjAbun_sub<-NSCLC_vsv_os_adjAbun_sub[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]

NSCLC_dsv_resp_adjAbun_sub<-NSCLC_dsv_resp_adjAbun_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]

NSCLC_dsv_os_adjAbun_sub<-NSCLC_dsv_os_adjAbun_select[,c("Pheno","X","fdr.p","HR","SE","HR_left","HR_right")]
NSCLC_dsv_os_adjAbun_sub$pvalue_meta<-NSCLC_dsv_os_adjAbun_sub$fdr.p
NSCLC_dsv_os_adjAbun_sub$het_p<-1
NSCLC_dsv_os_adjAbun_sub$effect_meta<-NSCLC_dsv_os_adjAbun_sub$HR
NSCLC_dsv_os_adjAbun_sub$cilb<-NSCLC_dsv_os_adjAbun_sub$HR_left
NSCLC_dsv_os_adjAbun_sub$ciub<-NSCLC_dsv_os_adjAbun_sub$HR_right
NSCLC_dsv_os_adjAbun_sub$se_meta<-NSCLC_dsv_os_adjAbun_sub$SE
NSCLC_dsv_os_adjAbun_sub<-NSCLC_dsv_os_adjAbun_sub[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]

NSCLC_vsv_resp_adjAbun_sub$SV<-c(paste("vsv:",NSCLC_vsv_resp_adjAbun_sub$X,sep = ""))
NSCLC_vsv_os_adjAbun_sub$SV<-c(paste("vsv:",NSCLC_vsv_os_adjAbun_sub$X,sep = ""))

NSCLC_dsv_resp_adjAbun_sub$SV<-c(paste("dsv:",NSCLC_dsv_resp_adjAbun_sub$X,sep = ""))
NSCLC_dsv_os_adjAbun_sub$SV<-c(paste("dsv:",NSCLC_dsv_os_adjAbun_sub$X,sep = ""))

sv_prog_meta_pvalue.sig.edge<-rbind(NSCLC_vsv_resp_adjAbun_sub,NSCLC_vsv_os_adjAbun_sub,
NSCLC_dsv_resp_adjAbun_sub,NSCLC_dsv_os_adjAbun_sub)

sv_prog_meta_pvalue.sig.edge$beta<-log(sv_prog_meta_pvalue.sig.edge$effect_meta)
sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge$Sign<-NA

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.1 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "*" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.05 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "☆" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.01 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "★" 

sv_prog_meta_pvalue.sig.edge<-unique(sv_prog_meta_pvalue.sig.edge)

sv_prog_adjAbun.r<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "beta")] %>% spread("Pheno", "beta") %>% data.frame(row.names = "SV")
sv_prog_adjAbun.sign<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "Sign")] %>% spread("Pheno", "Sign") %>% data.frame(row.names = "SV")

sv_prog_meta_pvalue.sig.edge<-subset(sv_prog_meta_pvalue.sig.edge,is.na(pvalue_meta)==F,drop=T)
# write.csv(sv_prog_meta_pvalue.sig.edge,"test.csv")


##########################################################################################################
###  A.muciniphila  Akkermansia muciniphila ATCC BAA-835:

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("muciniphila", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("muciniphila", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("muciniphila", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("muciniphila", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Akkermansia muciniphila ATCC BAA-835:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-2, 0, length.out=ceiling(i/2) + 1), 
              seq(2/i, 2, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_A.muciniphila_sv_lm_adjAbun.heatmap.tiff", width =2800, height =1800, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("NSCLC\nResponse","NSCLC\n OS"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =3,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()



##########################################################################################################
###  B.caccae

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("caccae", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("caccae", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("caccae", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("caccae", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Bacteroides caccae ATCC 43185:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-2, 0, length.out=ceiling(i/2) + 1), 
              seq(2/i, 2, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_B.caccae_sv_adjAbun.heatmap.tiff", width =2800, height =1800, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("NSCLC\nResponse","NSCLC\n OS"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =3,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
  )

dev.off()



################################################################################################
### RCC


NSCLC_vsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/NSCLC_vsv_resp_meta_pvalue.csv",header=T)
NSCLC_vsv_resp_meta_pvalue.anno.edge<-left_join(NSCLC_vsv_resp_meta_pvalue.edge, France_vsv_info, by = c("X" = "SV_Name"))

NSCLC_vsv_resp_adjAbun_sig_spe<-subset(NSCLC_vsv_resp_meta_pvalue.edge,as.numeric(pvalue_meta) <= 0.01 & het_p > 0.05 & p_count >=2,drop=T)$X

vsv_select<-data.frame(c(NSCLC_vsv_resp_adjAbun_sig_spe))
colnames(vsv_select)<-"X"

RCC_dsv_resp_meta_pvalue.edge<-read.csv("07.Microbial_GWAS/datasets/RCC_dsv_resp_pvalue_adjust.csv",header=T)
RCC_dsv_resp_meta_pvalue.anno.edge<-left_join(RCC_dsv_resp_meta_pvalue.edge, France_dsv_info, by = c("X" = "SV_Name"))

NSCLC_dsv_resp_adjAbun_sig_spe<-subset(NSCLC_dsv_resp_meta_pvalue.edge,as.numeric(pvalue_meta) < 0.01 & het_p > 0.05 & p_count >= 2,drop=T)$X
RCC_dsv_resp_adjAbun_sig_spe<-subset(RCC_dsv_resp_meta_pvalue.edge,as.numeric(fdr.p) <= 0.1,drop=T)$X

dsv_select<-data.frame(c(NSCLC_dsv_resp_adjAbun_sig_spe,RCC_dsv_resp_adjAbun_sig_spe))
colnames(dsv_select)<-"X"

NSCLC_vsv_resp_adjAbun_select<-merge(vsv_select,NSCLC_vsv_resp_meta_pvalue.anno.edge,by.x="X",by.y="X")

NSCLC_dsv_resp_adjAbun_select<-merge(dsv_select,NSCLC_dsv_resp_meta_pvalue.anno.edge,by.x="X",by.y="X")
RCC_dsv_resp_adjAbun_select<-merge(dsv_select,RCC_dsv_resp_meta_pvalue.anno.edge,by.x="X",by.y="X")

NSCLC_vsv_resp_adjAbun_select$Pheno<-"NSCLC_Response"
NSCLC_dsv_resp_adjAbun_select$Pheno<-"NSCLC_Response"
RCC_dsv_resp_adjAbun_select$Pheno<-"RCC_Response"

NSCLC_vsv_resp_adjAbun_sub<-NSCLC_vsv_resp_adjAbun_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]
NSCLC_dsv_resp_adjAbun_sub<-NSCLC_dsv_resp_adjAbun_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]

RCC_dsv_resp_adjAbun_sub<-RCC_dsv_resp_adjAbun_select[,c("Pheno","X","fdr.p","OR","SE","OR_left","OR_right")]
RCC_dsv_resp_adjAbun_sub$pvalue_meta<-RCC_dsv_resp_adjAbun_sub$fdr.p
RCC_dsv_resp_adjAbun_sub$het_p<-1
RCC_dsv_resp_adjAbun_sub$effect_meta<-RCC_dsv_resp_adjAbun_sub$OR
RCC_dsv_resp_adjAbun_sub$cilb<-RCC_dsv_resp_adjAbun_sub$OR_left
RCC_dsv_resp_adjAbun_sub$ciub<-RCC_dsv_resp_adjAbun_sub$OR_right
RCC_dsv_resp_adjAbun_sub$se_meta<-RCC_dsv_resp_adjAbun_sub$SE
RCC_dsv_resp_adjAbun_sub<-RCC_dsv_resp_adjAbun_sub[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","cilb","ciub")]

NSCLC_vsv_resp_adjAbun_sub$SV<-c(paste("vsv:",NSCLC_vsv_resp_adjAbun_sub$X,sep = ""))
NSCLC_dsv_resp_adjAbun_sub$SV<-c(paste("dsv:",NSCLC_dsv_resp_adjAbun_sub$X,sep = ""))

RCC_dsv_resp_adjAbun_sub$SV<-c(paste("dsv:",RCC_dsv_resp_adjAbun_sub$X,sep = ""))

sv_prog_meta_pvalue.sig.edge<-rbind(NSCLC_vsv_resp_adjAbun_sub,NSCLC_dsv_resp_adjAbun_sub,
RCC_dsv_resp_adjAbun_sub)

sv_prog_meta_pvalue.sig.edge$beta<-log(sv_prog_meta_pvalue.sig.edge$effect_meta)
sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge$Sign<-NA

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.1 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "*" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.05 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "☆" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.01 &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "★" 

sv_prog_meta_pvalue.sig.edge<-unique(sv_prog_meta_pvalue.sig.edge)

sv_prog_adjAbun.r<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "beta")] %>% spread("Pheno", "beta") %>% data.frame(row.names = "SV")
sv_prog_adjAbun.sign<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "Sign")] %>% spread("Pheno", "Sign") %>% data.frame(row.names = "SV")

sv_prog_meta_pvalue.sig.edge<-subset(sv_prog_meta_pvalue.sig.edge,is.na(pvalue_meta)==F,drop=T)
# write.csv(sv_prog_meta_pvalue.sig.edge,"test.csv")


##########################################################################################################
### P.distasonis

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("distasonis", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("distasonis", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("distasonis", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("distasonis", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Parabacteroides distasonis ATCC 8503:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-3, 0, length.out=ceiling(i/2) + 1), 
              seq(3/i, 3, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_P.distasonis_sv_adjAbun.heatmap.tiff", width =2800, height =1800, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("NSCLC\nResponse","RCC\nResponse"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =3,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
  )

dev.off()






