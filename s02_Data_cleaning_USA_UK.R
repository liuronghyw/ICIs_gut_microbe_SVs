### SV data processing
### 2023-09-17
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")

source("functions.R")
##install.packages("ggmosaic")
##library("plyr")
##library("dplyr")

#####################################################################
### Read SV files
dsgv<- read.delim("00.rawData/SV/USA_UK_dsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/USA_UK_vsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")

dsgv_anno<-read.delim("00.rawData/SV/s02.dSVs_USA_UK_anno.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
vsgv_anno<-read.delim("00.rawData/SV/s03.vSVs_USA_UK_anno.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)

### Read database files
taxa_length <- read.table("00.rawdata/database/Species_genome_size.tsv",
                        sep = "\t", header = T,check.names = F,stringsAsFactors = F)
taxonomy<- read.csv("00.rawdata/database/representatives.genomes.taxonomy.csv",
                   sep = ",", header = T,check.names = F,stringsAsFactors = F)
taxonomy[taxonomy == ""]<-"Unknown"
colnames(taxonomy)[1]<-'X'
ncbi<-read.csv("00.rawData/database/NCBI_accession.txt", sep = "\t",header = T)
tax_relationship<-read.csv("00.rawData/database/progenome1_species_relationship.tsv",sep = "\t",header = F)


############################################################################################
#### 1 Clean SV data
#### Get clean profiles
## overlap with clinical information and abdundance file 

clinical<-read.csv("00.rawData/clinical/clinical_whole_reads.csv",header=T)
rownames(clinical)<-clinical$id

SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)

clin_dsgv_inter <- intersect(rownames(clinical),rownames(dsgv))
clin_abun_inter <- intersect(clin_dsgv_inter,colnames(SV_abun_s)[-c(1,2)])

name1<-data.frame(clinical$id)
colnames(name1)<-"id"
name1$clin<-"T"

name2<-data.frame(rownames(dsgv))
colnames(name2)<-"id"
name2$dsgv<-"T"

name<-merge(name1,name2,by.x="id",by.y="id",all.x=T,all.y=T)
write.csv(name,"sample_clin_dSV_check.csv")

clinical_sub <- clinical[clin_abun_inter,]
write.table(clinical_sub, "01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv",sep = '\t')

dsgv<-dsgv[clin_abun_inter,]
vsgv<-vsgv[clin_abun_inter,]

# Change SV names
colnames(dsgv) <- changeSVname(colnames(dsgv))
colnames(vsgv) <- changeSVname(colnames(vsgv))

### the ratio of NAs
dsgv_non_NA_rate<-apply(dsgv[,-1], 2, myfun<-function(x){sum(is.na(x))/439})
vsgv_non_NA_rate<-apply(vsgv[,-1], 2, myfun<-function(x){sum(is.na(x))/439})

## Outputs
if(!dir.exists("01.cleanData")){dir.create("01.cleanData")}
if(!dir.exists("01.cleanData/SV_all")){dir.create("01.cleanData/SV_all")}

write.csv(dsgv_non_NA_rate,"01.cleanData/SV_info/USA_UK_dsgv_non_NA_rate.csv")
write.csv(vsgv_non_NA_rate,"01.cleanData/SV_info/USA_UK_vsgv_non_NA_rate.csv")


### 2 Get name conversion table
###  Name conversion
organism<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  .[!duplicated(.)]

Short_name<- organism %>% 
  str_replace_all('\\[','') %>%
  str_replace_all('\\]', '') %>%
  str_replace_all(' cf\\.','')

Short_name[grep(' sp\\.', organism, invert = F)] <- Short_name[grep(' sp\\.', organism, invert = F)] %>%
  str_replace_all('sp\\..*','sp')

Fst_letter<-Short_name[grep(' sp\\.', organism, invert = T)] %>%
  str_replace_all(' .*','') %>%
  str_sub(start = 1,end = 1)

Spe_name<-Short_name[grep(' sp\\.', organism, invert = T)] %>%
  str_extract_all(' .*') %>%
  str_replace_all('^ ', '') %>%
  str_replace_all(' .*', '')

Short_name[grep(' sp\\.', organism, invert = T)] <-paste(Fst_letter,'.', Spe_name, sep = '')
# write.csv(Short_name,"Short_name_test.csv")

taxa_name<-data.frame(NCBI_taxonomy_id = taxonomy$X[match(organism,taxonomy$organism)],
                      organism = as.character(organism), 
                      Short_name = as.character(Short_name), stringsAsFactors = F)

taxa_name$Short_name[match('bacterium LF-3',taxa_name$organism)]<-'bacterium LF-3'
taxa_name<-left_join(taxa_name, ncbi, by = "NCBI_taxonomy_id")

if(!dir.exists("01.cleanData/SV_info")){dir.create("01.cleanData/SV_info")}
#write.table(taxa_name, "01.cleanData/SV_info/Species_name.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

#write.csv(taxa_name, "01.cleanData/SV_info/Species_name.csv",row.names=F)
## fjy<-read.csv("01.cleanData/SV_info/article_fjy.csv",header=T)
## overlap<-merge(taxa_name,fjy,by.x="NCBI_taxonomy_id",by.y="NCBI_taxonomy_id")
## write.csv(overlap,"01.cleanData/SV_info/overlap_article_fjy.csv")

### 4.3 Get SV annotation tables
# SV annotation tables
dsgv_info_anno<-data.frame(dsgv_anno,
                           SV_ID=dsgv_anno$SV_id,
                           Taxonomy_Name = taxa_name$organism[match(str_replace_all(dsgv_anno$Taxonomy_id, '\\..*', ''),
                                                                    taxa_name$NCBI_taxonomy_id)],
                           SV_Name = changeSVname(dsgv_anno$SV_id),
                           Taxonomy_ID = dsgv_anno$Taxonomy_id,
                           SV_size = calcSVSize(dsgv_anno$SV_id))[,c(9,7,6,8,3,10,4,5)]

vsgv_info_anno<-data.frame(vsgv_anno,
                           SV_ID=vsgv_anno$SV_id,
                           Taxonomy_Name = taxa_name$organism[match(str_replace_all(vsgv_anno$Taxonomy_id, '\\..*', ''),
                                                                    taxa_name$NCBI_taxonomy_id)],
                           SV_Name = changeSVname(vsgv_anno$SV_id),
                           Taxonomy_ID = vsgv_anno$Taxonomy_id,
                           SV_size = calcSVSize(vsgv_anno$SV_id))[,c(9,7,6,8,3,10,4,5)]

write.table(dsgv_info_anno, "01.cleanData/SV_info/USA_UK_dsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)
write.table(vsgv_info_anno, "01.cleanData/SV_info/USA_UK_vsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)

write.csv(dsgv_info_anno, "01.cleanData/SV_info/USA_UK_dsgv_info_anno.csv",row.names = F)
write.csv(vsgv_info_anno, "01.cleanData/SV_info/USA_UK_vsgv_info_anno.csv",row.names = F)


###########################################################################################
### 4.4 Get species information table
###  Get SV number per species

species_dsgv_n<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  table(.) %>%
  as.data.frame(.)
colnames(species_dsgv_n)<-c("Species","Deletion SVs number")
species_vsgv_n<-str_replace_all(colnames(vsgv),"\\:\\d+_\\d+.*","") %>%
  table(.) %>%
  as.data.frame(.)
colnames(species_vsgv_n)<-c("Species","Variable SVs number")

species_sgv_n<-full_join(species_dsgv_n, species_vsgv_n, by = "Species")
species_sgv_n[is.na(species_sgv_n)]<-0

NCBI_taxonomy_id<-species_sgv_n$Species %>%
  match(.,taxonomy$organism) %>%
  taxonomy$X[.]
species_sgv_n<-data.frame(NCBI_taxonomy_id, species_sgv_n)

## Get sample size per species
dsgv_infor_sample_n<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  duplicated(.) %>%
  `!`%>%
  dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)
colnames(dsgv_infor_sample_n) <- "Sample_number"
rownames(dsgv_infor_sample_n) <- rownames(dsgv_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
dsgv_infor_sample_n<-data.frame(Species = rownames(dsgv_infor_sample_n),dsgv_infor_sample_n)

Taxonomy_name <- match(dsgv_infor_sample_n$Species,taxa_name$organism) %>%
  taxa_name$Short_name[.]
sample_n<-data.frame(Short_name=Taxonomy_name, dsgv_infor_sample_n)


### output different cohort
clinical <- read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv")
table(clinical$dataset)

PRJNA397906<-subset(clinical,dataset=="PRJNA397906",drop=T)$id
PRJNA762360<-subset(clinical,dataset=="PRJNA762360",drop=T)$id
PRJNA770295<-subset(clinical,dataset=="PRJNA770295",drop=T)$id
PRJEB43119<-subset(clinical,dataset=="PRJEB43119",drop=T)$id
PRJNA541981<-subset(clinical,dataset=="PRJNA541981",drop=T)$id

PRJNA397906_dsgv<- dsgv[rownames(dsgv) %in% PRJNA397906,] 
PRJNA762360_dsgv<- dsgv[rownames(dsgv) %in% PRJNA762360,] 
PRJNA770295_dsgv<- dsgv[rownames(dsgv) %in% PRJNA770295,] 
PRJEB43119_dsgv<- dsgv[rownames(dsgv) %in% PRJEB43119,]
PRJNA541981_dsgv<- dsgv[rownames(dsgv) %in% PRJNA541981,]

PRJNA397906_vsgv<- vsgv[rownames(vsgv) %in% PRJNA397906,] 
PRJNA762360_vsgv<- vsgv[rownames(vsgv) %in% PRJNA762360,] 
PRJNA770295_vsgv<- vsgv[rownames(vsgv) %in% PRJNA770295,] 
PRJEB43119_vsgv<- vsgv[rownames(vsgv) %in% PRJEB43119,]
PRJNA541981_vsgv<- vsgv[rownames(vsgv) %in% PRJNA541981,]


## PRJNA397906 sample size per species
PRJNA397906_infor_sample_n<-str_replace_all(colnames(PRJNA397906_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA397906_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA397906_infor_sample_n) <- "PRJNA397906"
rownames(PRJNA397906_infor_sample_n) <- rownames(PRJNA397906_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA397906_infor_sample_n<-data.frame(Species = rownames(PRJNA397906_infor_sample_n),PRJNA397906_infor_sample_n)


## PRJNA762360 sample size per species
PRJNA762360_infor_sample_n<-str_replace_all(colnames(PRJNA762360_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA762360_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA762360_infor_sample_n) <- "PRJNA762360"
rownames(PRJNA762360_infor_sample_n) <- rownames(PRJNA762360_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA762360_infor_sample_n<-data.frame(Species = rownames(PRJNA762360_infor_sample_n),PRJNA762360_infor_sample_n)


## PRJNA770295 sample size per species
PRJNA770295_infor_sample_n<-str_replace_all(colnames(PRJNA770295_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA770295_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA770295_infor_sample_n) <- "PRJNA770295"
rownames(PRJNA770295_infor_sample_n) <- rownames(PRJNA770295_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA770295_infor_sample_n<-data.frame(Species = rownames(PRJNA770295_infor_sample_n),PRJNA770295_infor_sample_n)


## PRJEB43119 sample size per species
PRJEB43119_infor_sample_n<-str_replace_all(colnames(PRJEB43119_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJEB43119_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJEB43119_infor_sample_n) <- "PRJEB43119"
rownames(PRJEB43119_infor_sample_n) <- rownames(PRJEB43119_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJEB43119_infor_sample_n<-data.frame(Species = rownames(PRJEB43119_infor_sample_n),PRJEB43119_infor_sample_n)


## PRJNA541981 sample size per species
PRJNA541981_infor_sample_n<-str_replace_all(colnames(PRJNA541981_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA541981_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA541981_infor_sample_n) <- "PRJNA541981"
rownames(PRJNA541981_infor_sample_n) <- rownames(PRJNA541981_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA541981_infor_sample_n<-data.frame(Species = rownames(PRJNA541981_infor_sample_n),PRJNA541981_infor_sample_n)


## merge sample size from different dataset 
info_sample<-merge(PRJNA397906_infor_sample_n,PRJNA762360_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,PRJNA770295_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,PRJEB43119_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,PRJNA541981_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,sample_n,by.x="Species",by.y="Species")

infor_sample_n <- info_sample[,c(1,7,8,seq(2,6))]

## Merge sample size and SV number information
species_sample_n<-dplyr::full_join(species_sgv_n,infor_sample_n, by = "Species")
taxa_length$Species<-str_replace_all(taxa_length$Species, '\\..*', '')
species_sample_n$NCBI_taxonomy_id<-as.character(species_sample_n$NCBI_taxonomy_id)
species_sample_n<-dplyr::left_join(species_sample_n, taxa_length, by = c("NCBI_taxonomy_id"="Species"))
species_sample_n<-data.frame(species_sample_n,
                             SVs.number = species_sample_n[,3]+species_sample_n[,4])

## Merge all information
Informative_species_information <- match(species_sample_n$NCBI_taxonomy_id, taxonomy$X)%>%
  taxonomy[.,] %>%
  cbind(.,species_sample_n)

taxa_name_short<-taxa_name[,c(2,4)]

info <- full_join(Informative_species_information[,-11],
                                             taxa_name_short,
                                             by = 'organism')[,c(1:9,13,22,11,12,21,seq(15,19),14,20)]

colnames(info)[c(1,10:21)]<-c("NCBI_taxonomy_id","Short_name","NCBI_bioproject_accession", 
"Deletion_SVs_number", "Variable_SVs_number","SVs_number","PRJNA397906_sample_number",
"PRJNA762360_sample_number","PRJNA770295_sample_number",
"PRJEB43119_sample_number","PRJNA541981_sample_number","Total_samples_number", "Length")

write.table(info, "01.cleanData/SV_info/USA_UK_Informative_species_information.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info, "01.cleanData/SV_info/USA_UK_Informative_species_information.csv", row.names = F)


##########################################################################
### 筛选species,remove the replicate
info<-subset(info,is.na(NCBI_taxonomy_id)==F,drop=T)
info_sub<-subset(info, organism!="Oscillibacter sp. KLE 1728" & is.na(Short_name)==F
     & organism!="Enterococcus faecium NRRL B-2354",drop=T)

#########################################################################
###  5 Clean species abundance data 

# Read species abundance files
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)
# Basic
all_basic <- read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv")

sv_spe_taxid<-info_sub$NCBI_taxonomy_id
sv_spe_taxid[sv_spe_taxid==245018]<-649756
sv_spe_taxid[sv_spe_taxid==1118061]<-2585118

# relative abundance of species detected with SVs
#SV_abun_s[,-c(1:2)]<-apply(SV_abun_s[,-c(1:2)], 2, myfun<-function(x){x/sum(x)})

SV_abun_s <- sv_spe_taxid %>% 
  match(., tax_relationship$V2) %>% 
  tax_relationship$V1[.] %>% 
  match(., SV_abun_s$taxonomy_id) %>% 
  SV_abun_s[., -c(1:2)] %>% 
  t %>%
  as.data.frame

colnames(SV_abun_s)<-info_sub$organism
#rownames(SV_abun_s)<-str_replace_all(rownames(SV_abun_s), ".S.bracken", "")

SV_abun_s<-SV_abun_s[match(rownames(all_basic),rownames(SV_abun_s)),]
rownames(SV_abun_s)<-rownames(all_basic)

spe_abun_mean<-apply(SV_abun_s,2,mean,na.rm=T)
write.csv(spe_abun_mean,"01.cleanData/spe_abun_mean.csv")

nr<-nrow(SV_abun_s)
nc<-ncol(SV_abun_s)

zero_ratio<-c()

for(j in 1:nc){
   count<-0
   for(i in 1:nr){
     if(SV_abun_s[i,j]==0){count=count+1}
     }
   zero_ratio[j]<-count/nr
  }

zero_ratio<-data.frame(zero_ratio)
zero_ratio$organism<-colnames(SV_abun_s)
#write.csv(zero_ratio,"01.cleanData/spe_abun_zero_ratio.csv")
zero_ratio<-subset(zero_ratio,zero_ratio<0.4,drop=T)

info_final<-merge(info_sub,zero_ratio,by.x="organism",by.y="organism",all.x=T)[,-22]
head(info_final)

write.table(info_final, "01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info_final, "01.cleanData/SV_info/USA_UK_Informative_species_information_final.csv", row.names = F)


############################### abundance final
# Read species abundance files
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)
# Basic
all_basic <- read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv")

sv_spe_taxid<-info_final$NCBI_taxonomy_id
sv_spe_taxid[sv_spe_taxid==245018]<-649756
sv_spe_taxid[sv_spe_taxid==1118061]<-2585118

# relative abundance of species detected with SVs
#SV_abun_s[,-c(1:2)]<-apply(SV_abun_s[,-c(1:2)], 2, myfun<-function(x){x/sum(x)})

SV_abun_s <- sv_spe_taxid %>% 
  match(., tax_relationship$V2) %>% 
  tax_relationship$V1[.] %>% 
  match(., SV_abun_s$taxonomy_id) %>% 
  SV_abun_s[., -c(1:2)] %>% 
  t %>%
  as.data.frame

colnames(SV_abun_s)<-info_final$organism
#rownames(SV_abun_s)<-str_replace_all(rownames(SV_abun_s), ".S.bracken", "")

SV_abun_s<-SV_abun_s[match(rownames(all_basic),rownames(SV_abun_s)),]
rownames(SV_abun_s)<-rownames(all_basic)

# outputs
if(!dir.exists("01.cleanData/mbio_all")){dir.create("01.cleanData/mbio_all")}
write.table(SV_abun_s, "01.cleanData/mbio_all/USA_UK_SV_species_abun.tsv",sep = '\t')

mean(rowSums(SV_abun_s,na.rm = T),na.rm = T)  # Mean of total relative abundance of species detected with SVs: 0.5003941
rowSums(SV_abun_s,na.rm = T)[rowSums(SV_abun_s,na.rm = T)!=0] %>% min # Minimum of total relative abundance of species detected with SVs: 0.02413
max(rowSums(SV_abun_s,na.rm = T),na.rm = T)   # Maximum of total relative abundance of species detected with SVs: 0.81953


###########################################################################
## 不同的数据集
PRJNA397906_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA397906,] 
PRJNA762360_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA762360,] 
PRJNA770295_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA770295,] 
PRJEB43119_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJEB43119,] 
PRJNA541981_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA541981,]

write.table(PRJNA397906_abun, "01.cleanData/mbio_all/PRJNA397906_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA762360_abun, "01.cleanData/mbio_all/PRJNA762360_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA770295_abun, "01.cleanData/mbio_all/PRJNA770295_SV_species_abun.tsv",sep = '\t')
write.table(PRJEB43119_abun, "01.cleanData/mbio_all/PRJEB43119_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA541981_abun, "01.cleanData/mbio_all/PRJNA541981_SV_species_abun.tsv",sep = '\t')

mean(rowSums(SV_abun_s,na.rm = T),na.rm = T)  # Mean of total relative abundance of species detected with SVs: 0.825154
rowSums(SV_abun_s,na.rm = T)[rowSums(SV_abun_s,na.rm = T)!=0] %>% min # Minimum of total relative abundance of species detected with SVs: 0.43729
max(rowSums(SV_abun_s,na.rm = T),na.rm = T)   # Maximum of total relative abundance of species detected with SVs: 0.94705


###################################################################################################
###  select dsgv and vsgv information
info_final<-read.delim("01.cleanData/SV_info/USA_UK_Informative_species_information_final.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
organism<-info_final$organism

select_column<-c()
for(i in 1:ncol(dsgv)){
  name<-str_replace_all(colnames(dsgv[i]),"\\:\\d+_\\d+.*","")
  if(name %in% organism){select_column<-c(select_column,i)}
  }
dsgv_sub<-dsgv[,select_column]

select_column<-c()
for(i in 1:ncol(vsgv)){
  name<-str_replace_all(colnames(vsgv[i]),"\\:\\d+_\\d+.*","")
  if(name %in% organism){select_column<-c(select_column,i)}
  }
vsgv_sub<-vsgv[,select_column]

### 当variation同时出现在dsv和vsv中，去掉vsv信息

dsv_id<-colnames(dsgv_sub)
vsv_id<-colnames(vsgv_sub)
overlap_id<- intersect(dsv_id,vsv_id)
vsv_id_left<-setdiff(vsv_id,overlap_id)

vsgv_sub<-vsgv_sub[,colnames(vsgv_sub) %in% vsv_id_left]
write.csv(overlap_id,"check_overlap_id.csv")

write.table(dsgv_sub,"01.cleanData/SV_all/USA_UK_dsgv_all.tsv",sep = '\t')
write.table(vsgv_sub,"01.cleanData/SV_all/USA_UK_vsgv_all.tsv",sep = '\t')
save(dsgv_sub, file = "01.cleanData/SV_all/USA_UK_dsgv.RData")
save(vsgv_sub, file = "01.cleanData/SV_all/USA_UK_vsgv.RData")


### output different cohort
clinical <- read.table("01.cleanData/phen/USA_UK_Clinical_basic_overlap.tsv")
table(clinical$dataset)

PRJNA397906<-subset(clinical,dataset=="PRJNA397906",drop=T)$id
PRJNA762360<-subset(clinical,dataset=="PRJNA762360",drop=T)$id
PRJNA770295<-subset(clinical,dataset=="PRJNA770295",drop=T)$id
PRJEB43119<-subset(clinical,dataset=="PRJEB43119",drop=T)$id
PRJNA541981<-subset(clinical,dataset=="PRJNA541981",drop=T)$id

PRJNA397906_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA397906,] 
PRJNA762360_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA762360,] 
PRJNA770295_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA770295,] 
PRJEB43119_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJEB43119,] 
PRJNA541981_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA541981,]

PRJNA397906_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA397906,] 
PRJNA762360_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA762360,] 
PRJNA770295_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA770295,] 
PRJEB43119_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJEB43119,] 
PRJNA541981_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA541981,] 

write.table(PRJNA397906_dsgv,"01.cleanData/SV_all/dsgv_PRJNA397906.tsv",sep = '\t')
write.table(PRJNA397906_vsgv,"01.cleanData/SV_all/vsgv_PRJNA397906.tsv",sep = '\t')
save(PRJNA397906_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA397906.RData")
save(PRJNA397906_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA397906.RData")

write.table(PRJNA762360_dsgv,"01.cleanData/SV_all/dsgv_PRJNA762360.tsv",sep = '\t')
write.table(PRJNA762360_vsgv,"01.cleanData/SV_all/vsgv_PRJNA762360.tsv",sep = '\t')
save(PRJNA762360_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA762360.RData")
save(PRJNA762360_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA762360.RData")

write.table(PRJNA770295_dsgv,"01.cleanData/SV_all/dsgv_PRJNA770295.tsv",sep = '\t')
write.table(PRJNA770295_vsgv,"01.cleanData/SV_all/vsgv_PRJNA770295.tsv",sep = '\t')
save(PRJNA770295_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA770295.RData")
save(PRJNA770295_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA770295.RData")

write.table(PRJEB43119_dsgv,"01.cleanData/SV_all/dsgv_PRJEB43119.tsv",sep = '\t')
write.table(PRJEB43119_vsgv,"01.cleanData/SV_all/vsgv_PRJEB43119.tsv",sep = '\t')
save(PRJEB43119_dsgv, file = "01.cleanData/SV_all/dsgv_PRJEB43119.RData")
save(PRJEB43119_vsgv, file = "01.cleanData/SV_all/vsgv_PRJEB43119.RData")

write.table(PRJNA541981_dsgv,"01.cleanData/SV_all/dsgv_PRJNA541981.tsv",sep = '\t')
write.table(PRJNA541981_vsgv,"01.cleanData/SV_all/vsgv_PRJNA541981.tsv",sep = '\t')
save(PRJNA541981_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA541981.RData")
save(PRJNA541981_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA541981.RData")


################################ 4.5 Get distance matrices
#### 4.5.1 All samples

## msv (vsv+dsv) distance
#sgv<-cbind(vsgv_sub, dsgv_sub)

#all_shared_vsv_dis<-shared_sv_dis_canberra(vsgv)
#all_shared_dsv_dis<-shared_sv_dis_jaccard(dsgv)

all_shared_vsv_dis<-shared_sv_dis_canberra(vsgv_sub)
all_shared_dsv_dis<-shared_sv_dis_jaccard(dsgv_sub)

all_shared_sv_dis<-average_matrix(all_shared_vsv_dis,all_shared_dsv_dis)
save(all_shared_sv_dis, file = "01.cleanData/SV_all/USA_UK_all_shared_sv_dis.RData")

################################################################
## ## SV distance matrices of all species (vsv+dsv) distance overall and for different datasets

dataset_msv_dist(vsgv_sub,dsgv_sub,info_final,
"01.cleanData/SV_all/distMat/USA_UK_all_msv_dist.RData",
"01.cleanData/SV_all/distMat/USA_UK_all_msv_dist_std.RData")

dataset_msv_dist(PRJNA397906_vsgv,PRJNA397906_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJNA397906_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJNA397906_msv_dist_std.RData")

dataset_msv_dist(PRJNA762360_vsgv,PRJNA762360_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJNA762360_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJNA762360_msv_dist_std.RData")

dataset_msv_dist(PRJNA770295_vsgv,PRJNA770295_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJNA770295_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJNA770295_msv_dist_std.RData")

dataset_msv_dist(PRJEB43119_vsgv,PRJEB43119_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJEB43119_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJEB43119_msv_dist_std.RData")

dataset_msv_dist(PRJNA541981_vsgv,PRJNA541981_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJNA541981_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJNA541981_msv_dist_std.RData")


##########################################################################
## difference between cohorts

### Read SV files
dsgv<- read.delim("00.rawData/SV/USA_UK_dsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/USA_UK_vsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")

clinical<-read.csv("00.rawData/clinical/clinical_whole_reads.csv",header=T)[,c("id","study")]
rownames(clinical)<-clinical$id

##########################################################
###  for dsgv (chisq test)

dsgv_ID<-colnames(dsgv)
dsgv$id<-rownames(dsgv)

ndsgv<-ncol(dsgv)-1
pvalue<-c()

#i=4

for(i in 1:ndsgv){
   dsgv_ex<-dsgv[,c(ncol(dsgv),i)]
   ID<-dsgv_ex[1,1]
   colnames(dsgv_ex)[2]<-"dsv"

   cor<-merge(clinical,dsgv_ex,by.clinical="id",by.exp="id")
   cor<-subset(cor,is.na(dsv)==F,drop=T)

   mytable <- xtabs(~study+dsv, data=cor)　　　　
   a<-chisq.test(mytable)　　
   pvalue[i]<-a$p.value
  }

pvalue<-as.matrix(pvalue,ndsgv,1)
rownames(pvalue)<-dsgv_ID

write.csv(pvalue,"diff_study_SV/USA_UK_dsv_chisq_test.csv")



##############################################################################
## for vsgv (k-s test)
library(rstatix)

vsgv_ID<-colnames(vsgv)
vsgv$id<-rownames(vsgv)

nvsgv<-ncol(vsgv)-1
pvalue<-c()

#i=4
for(i in 1:nvsgv){
   vsgv_ex<-vsgv[,c(ncol(vsgv),i)]
   ID<-vsgv_ex[1,1]
   colnames(vsgv_ex)[2]<-"vsv"

   cor<-merge(clinical,vsgv_ex,by.clinical="id",by.exp="id")
   cor<-subset(cor,is.na(vsv)==F,drop=T)
　　
   a<-kruskal_test(vsv~as.factor(study),data =cor)　　
   pvalue[i]<-a$p
  }

  pvalue<-as.matrix(pvalue,nvsgv,1)
  rownames(pvalue)<-vsgv_ID

write.csv(pvalue,"diff_study_SV/USA_UK_vsv_chisq_test.csv")


