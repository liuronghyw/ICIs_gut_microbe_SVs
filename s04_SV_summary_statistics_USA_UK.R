### SV information
### 2023-9-20
### LiuRong

library("ggplot2")
library("reshape2")	
library("RColorBrewer")

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")
source("functions.R")

## 1 Preparation
### 1.1 Import
### 1.2 Inputs
# SV
all_dsv <- read.table("01.cleanData/SV_all/USA_UK_dsgv_all.tsv",check.names = F)
all_vsv <- read.table("01.cleanData/SV_all/USA_UK_vsgv_all.tsv",check.names = F)

info    <- read.table("01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv",
                      sep = "\t", header = T, stringsAsFactors = F)

# abundance
all_abun_sv<-read.table("01.cleanData/mbio_all/USA_UK_SV_species_abun.tsv",check.names = F)

# Basic
all_basic <- read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv",header=T)

### correlation between read count and number of SVs in a sample 
samp_dsgv_num<-apply(all_dsv, 1, myfun<-function(x){sum(!is.na(x))})
samp_vsgv_num<-apply(all_vsv, 1, myfun<-function(x){sum(!is.na(x))})

samp_sv_num<-data.frame(samp_dsgv_num+samp_vsgv_num)
samp_sv_num$id<-rownames(samp_sv_num)
colnames(samp_sv_num)[1]<-"sv_number"

whole<-merge(all_basic,samp_sv_num,by.x="id",by.y="id")
cor.test(whole$sv_number,whole$read_count)

whole_sub<-subset(whole,read_count<10000000,drop=T)
cor.test(whole_sub$sv_number,whole_sub$read_count)

p_scatter<-ggplot(data=whole,aes(x=read_count/1000000, y=sv_number),color = "white", alpha = 0.2,size = 0.5)+
  #geom_smooth(data=whole ,aes(x=read_count/1000000, y=sv_number),method = "lm", color = "white", alpha = 0.2,size = 0.5)+
  geom_point(col="#b5182b")+
  xlab("Read count (million)")+
  ylab("Number of SVs")+xlim(0,50)+
  #scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

tiff(file = "pics/FS2_USA_UK_read_count_sv_number.tiff", width =1000, height =800, res =300) 
print(p_scatter)
dev.off()


## 2 SV summary
## 2.1 Species number
## Number of species with SVs

if(!dir.exists("04.SV_summary_statistics")){dir.create("04.SV_summary_statistics")}

dsgv_species_n<-sum(info$Deletion_SVs_number>0)
vsgv_species_n<-sum(info$Variable_SVs_number>0)
total_species_n<-sum(info$Variable_SVs_number>0 | info$Deletion_SVs_number>0)

species_n<-data.frame(item = rep("Informative species number", 3),
                      categories = c("SVs","Deletion SVs", "Variable SVs"),
                      value = c(total_species_n,dsgv_species_n, vsgv_species_n))

species_n$categories <- factor(species_n$categories, levels = species_n$categories)
species_n$categories <- factor(species_n$categories,levels(species_n$categories)[c(2,3,1)])

p_species_n<-ggplot(species_n, aes(x=categories, y=value,label = value))+
  geom_bar(aes(fill = categories),stat = 'identity')+
  geom_text(position = position_stack(vjust = 0.5), color = "white")+
  xlab(NULL)+
  ylab("Number of Informative species")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                    breaks=c("Deletion SVs", "Variable SVs", "SVs"),
                    labels=c("Deletion SVs", "Variable SVs", "SVs"),
                    values = mycolor3)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.key = element_rect(fill = NA), 
        legend.position = "none",
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

tiff(file = "pics/FS1_USA_UK_Informative_species_number.tiff", width =4000, height =3200, res =300) 
print(p_species_n)
dev.off()


### 2.1 SV number
#### 2.1.1 SV total number

dsgv_n<-sum(info$Deletion_SVs_number)
vsgv_n<-sum(info$Variable_SVs_number)
sgv_n<-dsgv_n+vsgv_n

sv_n<-data.frame(items = rep("SVs number", 2),
                      categories = c("Deletion SV", "Variable SV"),
                      value = c(dsgv_n, vsgv_n))

p_pie<-my_pie(sv_n, 'SVs number',mycol =rev(mycolor2_blue_yellow))
p_pie_nolabel<-my_pie_nolabel(sv_n, 'SVs number',mycol =rev(mycolor2_blue_yellow))

#pdf("pics/Total_SVs_number.pdf",height =5, width = 5)
#print(p_pie)
#dev.off()

tiff(file = "pics/F2_USA_UK_Total_SVs_number.tiff", width =2300, height =2300, res =300) 
print(p_pie)
dev.off()

tiff(file = "pics/F2_USA_UK_Total_SVs_number_nolabel.tiff", width =2300, height =2300, res =300) 
print(p_pie_nolabel)
dev.off()


#### 2.1.2 SV number per species

info_svs<-aggregate(info$SVs_number, by=list(Category= info$Short_name), FUN=sum)
info_svs<-data.frame(info_svs)
colnames(info_svs)<-c("Short_name","SVs_number")

species_sgv_n_order<- info_svs$Short_name[order(info_svs$SVs_number, decreasing = T)]

#species_sgv_n_order<- info$Short_name[order(info$SVs_number, decreasing = T)]
species_sgv_n_long<-gather(info[,c(10,12,13)], "type", "number", c(2:3))

p_sv_species<-ggplot(species_sgv_n_long, aes(x=Short_name, y=number, group = type, fill = type))+
  geom_bar(position = "stack",stat="identity")+
  xlab(NULL)+
  ylab("Number of SVs")+
  scale_x_discrete(limits = species_sgv_n_order)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                      breaks = c("Deletion_SVs_number", "Variable_SVs_number"),
                      labels = c("Deletion SVs              ", "Variable SVs"),
                      values = mycolor2_blue_yellow)+ggtitle("UK")+
  theme(plot.subtitle = element_text(vjust =1), 
        plot.caption = element_text(vjust = 1), 
        plot.margin=unit(c(1,1,1,1),"cm"),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 65, hjust = 1,size = 45,family ="sans",color="black"),
        axis.text.y = element_text(size = 45,family ="sans",color="black"),
        axis.title.y = element_text(size = 50,family ="sans",color="black"),
        plot.title = element_text(size = 60,family ="sans",color="black"),
        legend.position=c(0.2,0.9),
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size =50,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",size = 0.2),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

#pdf("pics/SVs_number_species.pdf", height = 8, width = 16)
#print(p_sv_species)
#dev.off()

tiff(file = "pics/F2_USA_UK_SVs_number_species.tiff", width =5000, height =2200, res =300) 
print(p_sv_species)
dev.off()


### 2.2 Sample size per species
infor_sample_n_order<- info$Short_name[order(info$Total_samples_number, decreasing = T)]
infor_sample_n <- info[,c(10,seq(15,19))]
infor_sample_n.long <- gather(infor_sample_n,'Cohort', 'Sample_size', c(2:6))

p_sample_n<-ggplot(infor_sample_n.long, aes(x=Short_name, y=Sample_size,group = Cohort))+
  geom_bar(aes(fill = Cohort),position = "stack",stat="identity")+
  xlab(NULL)+
  ylab("Number of samples")+
  scale_x_discrete(limits = infor_sample_n_order)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                      breaks = c("PRJNA397906_sample_number","PRJNA762360_sample_number","PRJNA541981_sample_number",
                                  "PRJNA770295_sample_number","PRJEB43119_sample_number"),
                      labels = c("FrankelAE_2017","McCullochJA_2022","PetersBA_2019",
                                "SpencerCN_2021","LeeKA_2022"),
                      values =c("#0072B5FF","#20854EFF","#E18727FF","#7876B1FF","#FFDC91FF"))+
 theme(plot.subtitle = element_text(vjust =1), 
        plot.caption = element_text(vjust = 1), 
        plot.margin=unit(c(0,1,1,1),"cm"),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 65, hjust = 1,size = 20,family ="sans",color="black"),
        axis.text.y = element_text(size = 15,family ="sans",color="black"),
        axis.title.y = element_text(size = 25,family ="sans",color="black"),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size =20,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",size = 0.2),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


#pdf("03.SV_summary_statistics/FS2_Samples_number_species.pdf", height = 6, width = 11) 
#print(p_sample_n)
#dev.off()

tiff(file = "pics/FS1_USA_UK_Samples_number_species.tiff", width =2600, height =1100, res =300) 
print(p_sample_n)
dev.off()


### 2.3 SV number factors
##Show relationships between SV numbers and genome size and sample size.
species_sample_n<-info[,c(1,10,15,11:15,18,19,20,21)]

length_sv_n_cor<-cor.test(species_sample_n$Length/1000, species_sample_n$SVs_number)
text_r<-paste("r=",round(length_sv_n_cor$estimate,digits = 3),"\np=",format(length_sv_n_cor$p.value,digits = 3),sep = "")

p_scatter_pie<-ggplot()+
  geom_smooth(data=species_sample_n ,aes(x=Length/10000, y=SVs_number),method = "lm", color = "white", alpha = 0.2,size = 0.5)+
  geom_scatterpie(data=species_sample_n ,aes(x=Length/10000, y=SVs_number,r=Total_samples_number/50),
                  cols=c("Deletion_SVs_number","Variable_SVs_number"),color=NA, alpha = 0.75)+
  coord_equal()+
  annotate("text", -Inf, Inf, label = c(text_r),hjust = -0.1, vjust = 1)+
  xlab("Genome size (10 kbp)")+
  ylab("Number of SVs")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name   = NULL,
                      breaks = c("Deletion_SVs_number", "Variable_SVs_number"),
                      labels = c("Deletion SVs              ", "Variable SVs"),
                      values = mycolor2_green_blue)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


#pdf("pics/Sample_n_SVs_n_scatter_pie.pdf", height = 5, width = 4)
#print(p_scatter_pie)
#dev.off()

tiff(file = "pics/FS1_USA_UK_Sample_n_SVs_n_scatter_pie.tiff", width =1200, height =800, res =300) 
print(p_scatter_pie)
dev.off()



### 3 Species abundance
sv_species_total_abun<-data.frame(id = rownames(all_abun_sv),
                                  Total_abundance=rowSums(all_abun_sv))

mean(sv_species_total_abun$Total_abundance,na.rm = T)
se(sv_species_total_abun$Total_abundance)
min(sv_species_total_abun$Total_abundance,na.rm = T)
max(sv_species_total_abun$Total_abundance,na.rm = T)

p_sv_abun_density<-ggplot(sv_species_total_abun,aes(x=Total_abundance), color = "#2EC4B6")+
  geom_density(color = "#b5182b", fill = "#34a186", alpha = 0.5)+
  geom_rug(color = "#2EC4B6",length = unit(0.05, "npc"))+
  geom_vline(xintercept  = mean(sv_species_total_abun$Total_abundance,na.rm = T), linetype = "dashed",color = "#ee6352")+
  ylab('Density')+
  xlab('Abundance')+ggtitle("USA/UK")+
  scale_y_continuous(expand = c(0,0))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=50,family ="sans"),
        axis.text.y = element_text(size=50,family ="sans"),
        axis.title.x = element_text(size=50,family ="sans"),
        axis.title.y = element_text(size=50,family ="sans"),
        plot.title = element_text(size=50,family ="sans"),
        legend.position = "none",
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

#pdf("pics/sv_abun_density.pdf", width = 3, height = 2)
#print(p_sv_abun_density)
#dev.off()

tiff(file = "pics/FS1_USA_UK_sv_abun_density.tiff", width =2200, height =1600, res =300) 
print(p_sv_abun_density)
dev.off()


## 4 SV correlation

vsv_corr <- rcorr(as.matrix(all_vsv))
vsv_corr.r<-vsv_corr$r
vsv_corr.p<-vsv_corr$P
vsv_corr.n<-vsv_corr$n
vsv_corr.r.edge<-melt(vsv_corr.r)
vsv_corr.p.edge<-melt(vsv_corr.p)
vsv_corr.n.edge<-melt(vsv_corr.n)

vsv_corr.edge<-cbind(vsv_corr.r.edge, vsv_corr.p.edge[,-c(1:2)], vsv_corr.n.edge[,-c(1:2)])
colnames(vsv_corr.edge)<-c("SV1", "SV2", "R", "P", "N")

save(vsv_corr.edge, file = "04.SV_summary_statistics/USA_UK_vsv_corr.edge.RData")

dsv_corr <- rcorr(as.matrix(all_dsv))
dsv_corr.r<-dsv_corr$r
dsv_corr.p<-dsv_corr$P
dsv_corr.n<-dsv_corr$n
dsv_corr.r.edge<-melt(dsv_corr.r)
dsv_corr.p.edge<-melt(dsv_corr.p)
dsv_corr.n.edge<-melt(dsv_corr.n)

dsv_corr.edge<-cbind(dsv_corr.r.edge, dsv_corr.p.edge[,-c(1:2)], dsv_corr.n.edge[,-c(1:2)])
colnames(dsv_corr.edge)<-c("SV1", "SV2", "R", "P", "N")

save(dsv_corr.edge, file = "04.SV_summary_statistics/USA_UK_dsv_corr.edge.RData")


all_vsv_M_plot.list<-list()
all_dsv_M_plot.list<-list()

for (i in c(1:nrow(info))){
#  i<-1
  
  file_name<-str_replace_all(info$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  
  
  vsv_i<-all_vsv[,grep(spe_name,colnames(all_vsv))]
  dsv_i<-all_dsv[,grep(spe_name,colnames(all_dsv))]
  
  if(info$organism[i]=="Bifidobacterium adolescentis"){
    p_vsv_M_i<-ggplot() + theme_void()
    
  }else{
      p_vsv_M_i<- ggcorr(vsv_i,method = c("pairwise"), alpha = 0)+
    ggtitle(info$Short_name[i])+
    theme(legend.position = "none",
          plot.title = element_text(size = 4))
    
  }

  p_dsv_M_i<- ggcorr(dsv_i,method = c("pairwise"), alpha = 0)+
    ggtitle(info$Short_name[i])+
    theme(legend.position = "none",
          plot.title = element_text(size = 4))
  
  all_vsv_M_plot.list[[i]] <- p_vsv_M_i
  all_dsv_M_plot.list[[i]] <- p_dsv_M_i
}

#pdf("pics/vsv_correlation.pdf")
#plot_grid(plotlist=all_vsv_M_plot.list)
#dev.off()

#pdf("pics/dsv_correlation.pdf")
#plot_grid(plotlist=all_dsv_M_plot.list)
#dev.off()

tiff(file = "pics/FS1_USA_UK_vsv_correlation.tiff", width =2000, height =2000, res =300) 
plot_grid(plotlist=all_vsv_M_plot.list)
dev.off()

tiff(file = "pics/FS1_USA_UK_dsv_correlation.tiff", width =2000, height =2000, res =300) 
plot_grid(plotlist=all_dsv_M_plot.list)
dev.off()




