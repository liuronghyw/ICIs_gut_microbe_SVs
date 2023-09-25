### SV data processing
### 2023-09-19
### LiuRong

setwd("D:/analysis_lab/pharmacomicrobiomics/1_ICB_metagenome_SV/pipeline/4_R")

source("functions.R")
##install.packages("ggmosaic")
##library("plyr")
##library("dplyr")

#####################################################################
### Read SV files
dsgv<- read.delim("00.rawData/SV/France_dsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/France_vsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")

dsgv_anno<-read.delim("00.rawData/SV/s02.dSVs_France_anno.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
vsgv_anno<-read.delim("00.rawData/SV/s03.vSVs_France_anno.tsv",
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
write.csv(name,"test/France_sample_clin_dSV_check.csv")

clinical_sub <- clinical[clin_abun_inter,]
write.table(clinical_sub, "01.cleanData/phen/France_Clinical_basic_overlap.tsv",sep = '\t')

dsgv<-dsgv[clin_abun_inter,]
vsgv<-vsgv[clin_abun_inter,]

# Change SV names
colnames(dsgv) <- changeSVname(colnames(dsgv))
colnames(vsgv) <- changeSVname(colnames(vsgv))

### the ratio of NAs
dsgv_non_NA_rate<-apply(dsgv[,-1], 2, myfun<-function(x){sum(is.na(x))/1011})
vsgv_non_NA_rate<-apply(vsgv[,-1], 2, myfun<-function(x){sum(is.na(x))/1011})

## Outputs
if(!dir.exists("01.cleanData")){dir.create("01.cleanData")}
if(!dir.exists("01.cleanData/SV_all")){dir.create("01.cleanData/SV_all")}

write.csv(dsgv_non_NA_rate,"01.cleanData/SV_info/France_dsgv_non_NA_rate.csv")
write.csv(vsgv_non_NA_rate,"01.cleanData/SV_info/France_vsgv_non_NA_rate.csv")


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

write.table(dsgv_info_anno, "01.cleanData/SV_info/France_dsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)
write.table(vsgv_info_anno, "01.cleanData/SV_info/France_vsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)

write.csv(dsgv_info_anno, "01.cleanData/SV_info/France_dsgv_info_anno.csv",row.names = F)
write.csv(vsgv_info_anno, "01.cleanData/SV_info/France_vsgv_info_anno.csv",row.names = F)


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
clinical <- read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")
table(clinical$dataset)

PRJEB22863<-subset(clinical,dataset=="PRJEB22863",drop=T)$id
PRJNA751792<-subset(clinical,dataset=="PRJNA751792",drop=T)$id

nsclc<-subset(clinical,cancer_type=="NSCLC",drop=T)$id
rcc<-subset(clinical,cancer_type=="RCC",drop=T)$id

PRJEB22863_dsgv<- dsgv[rownames(dsgv) %in% PRJEB22863,] 
PRJNA751792_dsgv<- dsgv[rownames(dsgv) %in% PRJNA751792,] 

PRJEB22863_vsgv<- vsgv[rownames(vsgv) %in% PRJEB22863,] 
PRJNA751792_vsgv<- vsgv[rownames(vsgv) %in% PRJNA751792,] 

## PRJNA22863 sample size per species
PRJEB22863_infor_sample_n<-str_replace_all(colnames(PRJEB22863_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJEB22863_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJEB22863_infor_sample_n) <- "PRJEB22863"
rownames(PRJEB22863_infor_sample_n) <- rownames(PRJEB22863_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJEB22863_infor_sample_n<-data.frame(Species = rownames(PRJEB22863_infor_sample_n),PRJEB22863_infor_sample_n)


## PRJNA751792 sample size per species
PRJNA751792_infor_sample_n<-str_replace_all(colnames(PRJNA751792_dsgv),"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  PRJNA751792_dsgv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(PRJNA751792_infor_sample_n) <- "PRJNA751792"
rownames(PRJNA751792_infor_sample_n) <- rownames(PRJNA751792_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
PRJNA751792_infor_sample_n<-data.frame(Species = rownames(PRJNA751792_infor_sample_n),PRJNA751792_infor_sample_n)


## merge sample size from different dataset 

info_sample<-merge(PRJEB22863_infor_sample_n,PRJNA751792_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,sample_n,by.x="Species",by.y="Species")

infor_sample_n <- info_sample[,c(1,4,5,2,3)]

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
                                             by = 'organism')[,c(1:9,13,19,11,12,18,15,16,14,17)]

colnames(info)[c(1,10:18)]<-c("NCBI_taxonomy_id","Short_name","NCBI_bioproject_accession", 
"Deletion_SVs_number", "Variable_SVs_number","SVs_number","PRJEB22863_sample_number",
"PRJNA751792_sample_number","Total_samples_number", "Length")

write.table(info, "01.cleanData/SV_info/France_Informative_species_information.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info, "01.cleanData/SV_info/France_Informative_species_information.csv", row.names = F)

##########################################################################
### 筛选species,remove the replicate

info_sub<-subset(info, organism!="Oscillibacter sp. KLE 1728"
   & organism!="Phascolarctobacterium sp. CAG:266",drop=T)


#########################################################################
###  5 Clean species abundance data 

# Read species abundance files
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)
# Basic
all_basic <- read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")

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

info_final<-merge(info_sub,zero_ratio,by.x="organism",by.y="organism",all.x=T)[,-24]
head(info_final)

write.table(info_final, "01.cleanData/SV_info/France_Informative_species_information_final.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info_final, "01.cleanData/SV_info/France_Informative_species_information_final.csv", row.names = F)


############################### abundance final

# Read species abundance files
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)
# Basic
all_basic <- read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")

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
write.table(SV_abun_s, "01.cleanData/mbio_all/France_SV_species_abun.tsv",sep = '\t')

mean(rowSums(SV_abun_s,na.rm = T),na.rm = T)  # Mean of total relative abundance of species detected with SVs: 0.825154
rowSums(SV_abun_s,na.rm = T)[rowSums(SV_abun_s,na.rm = T)!=0] %>% min # Minimum of total relative abundance of species detected with SVs: 0.43729
max(rowSums(SV_abun_s,na.rm = T),na.rm = T)   # Maximum of total relative abundance of species detected with SVs: 0.94705


###########################################################################
## 不同的数据集

PRJEB22863_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJEB22863,] 
PRJNA751792_abun<- SV_abun_s[rownames(SV_abun_s) %in% PRJNA751792,] 

write.table(PRJEB22863_abun, "01.cleanData/mbio_all/PRJEB22863_SV_species_abun.tsv",sep = '\t')
write.table(PRJNA751792_abun, "01.cleanData/mbio_all/PRJNA751792_SV_species_abun.tsv",sep = '\t')

mean(rowSums(SV_abun_s,na.rm = T),na.rm = T)  # Mean of total relative abundance of species detected with SVs: 0.825154
rowSums(SV_abun_s,na.rm = T)[rowSums(SV_abun_s,na.rm = T)!=0] %>% min # Minimum of total relative abundance of species detected with SVs: 0.43729
max(rowSums(SV_abun_s,na.rm = T),na.rm = T)   # Maximum of total relative abundance of species detected with SVs: 0.94705


###################################################################################################
###  select dsgv and vsgv information

info_final<-read.delim("01.cleanData/SV_info/France_Informative_species_information_final.tsv",
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

write.table(dsgv_sub,"01.cleanData/SV_all/France_dsgv_all.tsv",sep = '\t')
write.table(vsgv_sub,"01.cleanData/SV_all/France_vsgv_all.tsv",sep = '\t')
save(dsgv_sub, file = "01.cleanData/SV_all/France_dsgv.RData")
save(vsgv_sub, file = "01.cleanData/SV_all/France_vsgv.RData")


### output different cohort
clinical <- read.table("01.cleanData/phen/France_Clinical_basic_overlap.tsv")
table(clinical$dataset)

PRJEB22863<-subset(clinical,dataset=="PRJEB22863",drop=T)$id
PRJEB22863_NSCLC<-subset(clinical,dataset=="PRJEB22863" & cancer_type=="NSCLC",drop=T)$id
PRJEB22863_RCC<-subset(clinical,dataset=="PRJEB22863" & cancer_type=="RCC",drop=T)$id
PRJNA751792<-subset(clinical,dataset=="PRJNA751792",drop=T)$id

PRJEB22863_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJEB22863,] 
PRJEB22863_NSCLC_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJEB22863_NSCLC,] 
PRJEB22863_RCC_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJEB22863_RCC,] 
PRJNA751792_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% PRJNA751792,] 

PRJEB22863_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJEB22863,] 
PRJEB22863_NSCLC_vsgv<- dsgv_sub[rownames(vsgv_sub) %in% PRJEB22863_NSCLC,] 
PRJEB22863_RCC_vsgv<- dsgv_sub[rownames(vsgv_sub) %in% PRJEB22863_RCC,] 
PRJNA751792_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% PRJNA751792,] 

nsclc_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% nsclc,] 
rcc_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% rcc,] 

nsclc_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% nsclc,] 
rcc_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% rcc,] 

write.table(PRJEB22863_dsgv,"01.cleanData/SV_all/dsgv_PRJEB22863.tsv",sep = '\t')
write.table(PRJEB22863_vsgv,"01.cleanData/SV_all/vsgv_PRJEB22863.tsv",sep = '\t')
save(PRJEB22863_dsgv, file = "01.cleanData/SV_all/dsgv_PRJEB22863.RData")
save(PRJEB22863_vsgv, file = "01.cleanData/SV_all/vsgv_PRJEB22863.RData")

write.table(PRJEB22863_NSCLC_dsgv,"01.cleanData/SV_all/dsgv_NSCLC_PRJEB22863.tsv",sep = '\t')
write.table(PRJEB22863_NSCLC_vsgv,"01.cleanData/SV_all/vsgv_NSCLC_PRJEB22863.tsv",sep = '\t')
save(PRJEB22863_NSCLC_dsgv, file = "01.cleanData/SV_all/dsgv_NSCLC_PRJEB22863.RData")
save(PRJEB22863_NSCLC_vsgv, file = "01.cleanData/SV_all/vsgv_NSCLC_PRJEB22863.RData")

write.table(PRJEB22863_RCC_dsgv,"01.cleanData/SV_all/dsgv_RCC_PRJEB22863.tsv",sep = '\t')
write.table(PRJEB22863_RCC_vsgv,"01.cleanData/SV_all/vsgv_RCC_PRJEB22863.tsv",sep = '\t')
save(PRJEB22863_RCC_dsgv, file = "01.cleanData/SV_all/dsgv_RCC_PRJEB22863.RData")
save(PRJEB22863_RCC_vsgv, file = "01.cleanData/SV_all/vsgv_RCC_PRJEB22863.RData")

write.table(PRJNA751792_dsgv,"01.cleanData/SV_all/dsgv_PRJNA751792.tsv",sep = '\t')
write.table(PRJNA751792_vsgv,"01.cleanData/SV_all/vsgv_PRJNA751792.tsv",sep = '\t')
save(PRJNA751792_dsgv, file = "01.cleanData/SV_all/dsgv_PRJNA751792.RData")
save(PRJNA751792_vsgv, file = "01.cleanData/SV_all/vsgv_PRJNA751792.RData")


################################ 4.5 Get distance matrices
#### 4.5.1 All samples

## add_matrix
## example
#m1 <- matrix(c(1:9),nrow=3,ncol=3,dimnames=list(c("c1","c2","c3"),c("c1","c2","c3")))
#m2 <- matrix(c(5:13),nrow=3,ncol=3,dimnames=list(c("c1","c4","c3"),c("c1","c4","c3")))
#merge_mat<-average_matrix(m1,m2)

## msv (vsv+dsv) distance
#sgv<-cbind(vsgv_sub, dsgv_sub)

all_shared_vsv_dis<-shared_sv_dis_canberra(vsgv_sub)
all_shared_dsv_dis<-shared_sv_dis_jaccard(dsgv_sub)

all_shared_sv_dis<-average_matrix(all_shared_vsv_dis,all_shared_dsv_dis)
save(all_shared_sv_dis, file = "01.cleanData/SV_all/France_all_shared_sv_dis.RData")


################################################################
## ## SV distance matrices of all species (vsv+dsv) distance overall and for different datasets

dataset_msv_dist(vsgv_sub,dsgv_sub,info_final,
"01.cleanData/SV_all/distMat/France_all_msv_dist.RData",
"01.cleanData/SV_all/distMat/France_all_msv_dist_std.RData")

dataset_msv_dist(PRJEB22863_vsgv,PRJEB22863_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJEB22863_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJEB22863_msv_dist_std.RData")

dataset_msv_dist(PRJEB22863_NSCLC_vsgv,PRJEB22863_NSCLC_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJEB22863_NSCLC_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJEB22863_NSCLC_msv_dist_std.RData")

dataset_msv_dist(PRJEB22863_RCC_vsgv,PRJEB22863_RCC_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJEB22863_RCC_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJEB22863_RCC_msv_dist_std.RData")

dataset_msv_dist(PRJNA751792_vsgv,PRJNA751792_dsgv,info_final,
"01.cleanData/SV_all/distMat/PRJNA751792_msv_dist.RData",
"01.cleanData/SV_all/distMat/PRJNA751792_msv_dist_std.RData")

dataset_msv_dist(nsclc_vsgv,nsclc_dsgv,info_final,
"01.cleanData/SV_all/distMat/nsclc_msv_dist.RData",
"01.cleanData/SV_all/distMat/nsclc_msv_dist_std.RData")

dataset_msv_dist(rcc_vsgv,rcc_dsgv,info_final,
"01.cleanData/SV_all/distMat/rcc_msv_dist.RData",
"01.cleanData/SV_all/distMat/rcc_msv_dist_std.RData")



##########################################################################
## difference between cohorts

### Read SV files
dsgv<- read.delim("00.rawData/SV/France_dsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/France_vsgv.csv",
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

write.csv(pvalue,"diff_study_SV/France_dsv_chisq_test.csv")



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

write.csv(pvalue,"diff_study_SV/France_vsv_chisq_test.csv")


