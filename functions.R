## color panels
mycolor2_blue_yellow <- c("#0072B5FF","#FFDC91FF")
mycolor2_green_blue  <- c("#20854EFF","#7876B1FF")
mycolor2_blue_red    <- c("#6F99ADFF","#BC3C29FF")
mycolor2_red_blue    <- c("#BC3C29FF","#6F99ADFF")
mycolor2_yellow_green   <- c("#E18727FF","#20854EFF")
mycolor2_pfs<-c("#EE4C97FF","#0072B5FF")

mycolor3 <- c("#7876B1FF","#6F99ADFF","#FFDC91FF")   
mycolor3_red_blue_yellow <- c("#BC3C29FF","#0072B5FF","#E18727FF")
mycolor3_green_red_yellow<- c("#20854EFF","#BC3C29FF","#E18727FF")

mycolor6<- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")
mycolor7<- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF")
mycolor8<- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")


#install.packages("ggmosaic")
#install.packages("tidyverse")
#install.packages("GGally")
#install.packages("ggthemes")
#install.packages("gplots")
#install.packages("VennDiagram")
#install.packages("circlize")
#install.packages("viridis")
#install.packages("wesanderson")
#install.packages("ggalluvial")
#install.packages(c("Hmisc","vegan","mediation"))
#install.packages("rlang",version="1.0.3")
#install.packages("glue",version="1.3.2")
#install.packages("tibble",version="3.0.0")
#install.packages("dplyr")
#install.packages("showtext")
#install.packages("Cairo")
#install.packages("scatterpie")
#install.packages("NbClust")
#install.packages("fpc")
#install.packages("tsne")
#install.packages("meta")

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.12")

### library(BiocManager)
### BiocManager::install("microbiome")


## libraty packages
library(tidyverse)
library(ggpubr)
library(GGally)
library(ggmosaic)
library(cowplot)
library(ggsci)
library(ggthemes)
library(gplots)
library(VennDiagram)
library(circlize)
library(viridis)
library(wesanderson)
library(ggalluvial)
library(meta)

library(reshape2)
library(Hmisc)
library(vegan)
library(microbiome)
library(mediation)
library(metafor)

## Convert SV name
changeSVname<-function(SVrawid){
  testname     <- SVrawid
  species_name <- as.character(taxonomy$organism[match(str_replace_all(testname, "\\..*",""), taxonomy$X)])
  region       <- str_replace_all(testname, ".*\\:","") 
  region_list  <- str_split(region,";") 
  
  region_suf   <- NULL
  i<-1
  for (i in c(1:length(region_list))){
    if(length(region_list[[i]]) <= 2){
      region_suf[i] <- paste(":",region[i],sep = "")
    }else{
      region_suf[i] <- paste(":",region_list[[i]][1]," and ",length(region_list[[i]])-1," segments", sep = "")
    }
    i <- i+1
  }
  paste(species_name,region_suf,sep = "")
}


## Calculate SV size
calcSVSize <- function(SVrawid){
  testname <- SVrawid
  region   <- str_replace_all(testname, ".*\\:","") 
  region_list <- str_split(region,";") 
  
  sv_size  <- NULL
  for (i in c(1:length(region_list))) {
    frag_df    <- str_split(region_list[[i]], "_") %>% unlist %>% as.numeric  %>% matrix(byrow = T, ncol = 2) %>% as.data.frame 
    sv_size[i] <- sum(as.numeric(frag_df$V2)-as.numeric(frag_df$V1))
  }
  sv_size
}



shared_sv_dis_canberra<-function(inDf){
  #inDf<-all_sv
  inList <- split(as.matrix(inDf), seq(nrow(inDf)))
  
  shared_n_func<-function(inVec,Vec_i){
    #inVec<-inList[[3]]
    Vec_i<-Vec_i
    insvdf<-data.frame(inVec,Vec_i) %>% na.omit
    sv_dis<- vegdist(t(insvdf), method = "canberra")
    #length(na.omit(Vec_i+inVec))
    return(sv_dis)
  }
  
  marker_n_mat<-matrix(NA, nrow = nrow(inDf), ncol = nrow(inDf))
  for (i in 1:nrow(inDf)) { #nrow(inDf)
    #i<-2
    Vec_i<-inList[[i]]
    shared_n_i<-sapply(inList, shared_n_func,Vec_i = Vec_i)
    marker_n_mat[i,]<-shared_n_i
  }
  rownames(marker_n_mat)<-rownames(inDf)
  colnames(marker_n_mat)<-rownames(inDf)
  
  return(marker_n_mat)
}


shared_sv_dis_jaccard<-function(inDf){
  #inDf<-all_sv
  inList <- split(as.matrix(inDf), seq(nrow(inDf)))
  
  shared_n_func<-function(inVec,Vec_i){
    #inVec<-inList[[3]]
    Vec_i<-Vec_i
    insvdf<-data.frame(inVec,Vec_i) %>% na.omit
    sv_dis<- vegdist(t(insvdf), method = "jaccard")
    #length(na.omit(Vec_i+inVec))
    return(sv_dis)
  }
  
  marker_n_mat<-matrix(NA, nrow = nrow(inDf), ncol = nrow(inDf))
  for (i in 1:nrow(inDf)) { #nrow(inDf)
    #i<-2
    Vec_i<-inList[[i]]
    shared_n_i<-sapply(inList, shared_n_func,Vec_i = Vec_i)
    marker_n_mat[i,]<-shared_n_i
  }
  rownames(marker_n_mat)<-rownames(inDf)
  colnames(marker_n_mat)<-rownames(inDf)
  
  return(marker_n_mat)
}







## Calculate SE
se <- function(x) {x<-na.omit(x);sqrt(var(x)/length(x))}


## Summary statistics
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N            = length2(xx[[col]], na.rm=na.rm),
                     nonZero_N    = sum(xx[[col]]!=0,na.rm = na.rm),
                     nonZero_rate = sum(xx[[col]]!=0,na.rm = na.rm)/length2(xx[[col]], na.rm=na.rm),
                     mean         = mean   (xx[[col]], na.rm=na.rm),
                     sd           = sd     (xx[[col]], na.rm=na.rm),
                     mindata      = min(xx[[col]], na.rm=na.rm),
                     maxdata      = max(xx[[col]], na.rm=na.rm),
                     quant1_  = quantile(xx[[col]], na.rm=na.rm)[2] ,
                     quant2_  = quantile(xx[[col]], na.rm=na.rm)[4] ,
                     quantmid = quantile(xx[[col]], na.rm=na.rm)[3] ,
                     mea      = mean   (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  
  ciMult   <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  datac$xd <- datac$mea-datac$se
  datac$xu<- datac$mea+datac$se
  return(datac)
}

## Remove samples with NA in distance matrix
dist_rmna<-function(inDist){
  while(sort(colSums(is.na(inDist)),decreasing = T)[1] > 0){
    rmid<-names(sort(colSums(is.na(inDist)),decreasing = T)[1])
    inDist<-inDist[-match(rmid, rownames(inDist)),-match(rmid, colnames(inDist))]
  }
  return(inDist)
}


## pie chart
my_pie<-function(in1,item,mycol = c("#F24D4D","#4D7EAE")){
  require(tibble)
  require(showtext)
  require(Cairo)
  require(ggsci)
  require(tibble)
  require(scales)
  require(ggrepel)
  require(forcats)
  require(scatterpie)
  
  showtext_auto()
  
  in1.tibble            <- as_tibble(subset(in1, in1$items%in%item))
  in1.tibble$categories <- fct_reorder(in1.tibble$categories, in1.tibble$value)
  in1.tibble            <- in1.tibble[order(in1.tibble$value, decreasing = TRUE), ]
  piepercent            <- round(100*in1.tibble$value/sum(in1.tibble$value), 2)
  my_labels             <- tibble(x.breaks = seq(1.1, 1.1, length.out = length(piepercent)),
                                  y.breaks = cumsum(in1.tibble$value) - in1.tibble$value/2,
                                  labels = paste(in1.tibble$categories, "\n",in1.tibble$value,", ",piepercent, "%", sep = ""),
                                  categories = in1.tibble$categories)
  
  pdfName    <- paste(item, ".pie.pdf", sep = "")
  ggplot(in1.tibble, aes(x = 1, y = value, fill = categories)) +
    ggtitle(paste(item)) +
    geom_bar(stat="identity", color='white') + 
    coord_polar(theta='y') + 
    theme(legend.position = "None",
          axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels.
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    #    scale_fill_brewer(palette = "Set3", direction = -1)+
    scale_fill_manual(values=mycol) + # set fill color manually
    geom_text_repel(data = my_labels, 
                    aes(x = x.breaks, y = y.breaks, label = labels),
                    size = 7,        # label text size
                    show.legend = FALSE,
                    inherit.aes = FALSE,
                    arrow = arrow(length=unit(0.01, "npc")),
                    force = 1,colour = "white",
                    segment.color = 'grey50'   ## grey50
    )
}




## pie chart
my_pie_nolabel<-function(in1,item,mycol = c("#F24D4D","#4D7EAE")){
  require(tibble)
  require(showtext)
  require(Cairo)
  require(ggsci)
  require(tibble)
  require(scales)
  require(ggrepel)
  require(forcats)
  require(scatterpie)
  
  showtext_auto()
  
  in1.tibble            <- as_tibble(subset(in1, in1$items%in%item))
  in1.tibble$categories <- fct_reorder(in1.tibble$categories, in1.tibble$value)
  in1.tibble            <- in1.tibble[order(in1.tibble$value, decreasing = TRUE), ]
  piepercent            <- round(100*in1.tibble$value/sum(in1.tibble$value), 2)
  my_labels             <- tibble(x.breaks = seq(1.1, 1.1, length.out = length(piepercent)),
                                  y.breaks = cumsum(in1.tibble$value) - in1.tibble$value/2,
                                  labels = paste(in1.tibble$categories, "\n",in1.tibble$value,", ",piepercent, "%", sep = ""),
                                  categories = in1.tibble$categories)
  
  pdfName    <- paste(item, ".pie.pdf", sep = "")
  ggplot(in1.tibble, aes(x = 1, y = value, fill = categories)) +
    ggtitle(paste(item)) +
    geom_bar(stat="identity", color='white') + 
    coord_polar(theta='y') + 
    theme(legend.position = "None",
          axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels.
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    #    scale_fill_brewer(palette = "Set3", direction = -1)+
    scale_fill_manual(values=mycol) #+ # set fill color manually
    #geom_text_repel(data = my_labels, 
    #                aes(x = x.breaks, y = y.breaks, label = labels),
    #                size = 16,        # label text size
    #                show.legend = FALSE,
    #                inherit.aes = FALSE,
    #                arrow = arrow(length=unit(0.01, "npc")),
    #                force = 1,colour = "black",
    #                segment.color = 'grey50'   ## grey50
    #)
}






## Batch adonis
my_adonis_terms<-function(inDist, inDf, covDf=NULL, covar=NULL){
  require(vegan)
  ## test
  #inDist <- distList[[i]]
  #inDf   <- inDf
  #covDf  <- covDf
  #covar  <- inCovar
  ## test
  
  covDf<- covDf[,covar]
  
  inDist_rmna <- dist_rmna(inDist)
  covDf_rmna   <- na.omit(covDf)
  
  my_adonis<-function(inMat, inVec, covdf = NULL){
    ## test
    #inMat<-inDist_rmna
    #inVec<-inDf[,1]
    #covdf<-covDf_rmna
    ## test
    inMat_covdf_inter <- intersect(rownames(inMat),rownames(covdf))
    inMat_covdf_inVec_inter <- intersect(rownames(inDf)[!is.na(inVec)], inMat_covdf_inter)
    
    inMat_rmna <- inMat[match(inMat_covdf_inVec_inter,rownames(inMat)),
                        match(inMat_covdf_inVec_inter,colnames(inMat))]
    inVec_rmna <- inVec[match(inMat_covdf_inVec_inter, rownames(inDf))]
    covdf_rmna <- covdf[match(inMat_covdf_inVec_inter, rownames(covdf)),]
    
    in_cov <- data.frame(inVec_rmna,covdf_rmna)
    
    sample_size<-nrow(in_cov)

    if(length(table(inVec_rmna))<2){
      return(list(NA,NA, sample_size))
       }
    else{
      adonis_res <- adonis(as.dist(inMat_rmna)~.,data = in_cov)
      return(list(adonis_res$aov.tab[1,5],adonis_res$aov.tab[1,6], sample_size))
     }
  }
  
  adonis_res_batch <- apply(inDf, 2, my_adonis,inMat = inDist_rmna, covdf = covDf_rmna)
  adonis_res_batch.unlist <- matrix(unlist(adonis_res_batch), ncol = 3, byrow = T)
  colnames(adonis_res_batch.unlist)<-c("R2", "P", "N")
  
  return(adonis_res_batch.unlist)
}



my_adonis_terms_noadjAbun<-function(distList, inDf, covDf, covar,info){
  #distList<-PRJEB22863_msv_dist_std
  #inDf<-covar_PRJEB22863
  #covDf<-PRJEB22863_covar
  #covar<-covar2
  #info<-France_info
  
  res_table<-NULL
  for (i in 1:length(distList)) { # 
    #i<-8
    cat(paste(i,"\n")) 
    
    n<-names(distList)[i] %>% 
      str_replace_all("msv_","") %>% 
      match(., info$organism)
    
   inCovar<-covar
   #inCovar<-c(covar, info$organism[n])

    if(is.na(distList[[i]])==F){
      if (nrow(distList[[i]])>10 & ncol(distList[[i]])>10){
         adonis_res<-my_adonis_terms(distList[[i]], inDf, covDf = covDf, covar = inCovar)
          adonis_res_df<-data.frame(Species = rep(info$Short_name[n], dim(adonis_res)[1]),
                              Prog = colnames(inDf),
                              as.data.frame(adonis_res))
    
        res_table<-rbind(res_table, adonis_res_df)
      }
      else {
        adonis_res_df<-data.frame(Species = rep(info$Short_name[n], ncol(inDf)),
                              Prog = colnames(inDf),
                              R2=NA,P=NA,N=NA)
        res_table<-rbind(res_table, adonis_res_df)
       }
     }
    else {
        adonis_res_df<-data.frame(Species = rep(info$Short_name[n], ncol(inDf)),
                              Prog = colnames(inDf),
                              R2=NA,P=NA,N=NA)
        res_table<-rbind(res_table, adonis_res_df)
     }
  }
  
  res_table$Species<-as.character(res_table$Species)
  res_table$Prog<-as.character(res_table$Prog)
  
  res_table$fdr<-p.adjust(res_table$P,method = 'fdr')
  res_table$bonferroni<-p.adjust(res_table$P,method = 'bonferroni')
  
  # R2
  r2_mat<-matrix(data  = res_table$R2,
                 nrow  = length(unique(res_table$Species)),
                 ncol  = length(unique(res_table$Prog)),
                 byrow = T)
  
  rownames(r2_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(r2_mat)<-colnames(inDf)
  
  # P
  p_mat<-matrix(data  = res_table$P,
                nrow  = length(unique(res_table$Species)),
                ncol  = length(unique(res_table$Prog)),
                byrow = T)
  
  rownames(p_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(p_mat)<-colnames(inDf)
  
  # fdr
  fdr_mat<-matrix(data  = res_table$fdr,
                  nrow  = length(unique(res_table$Species)),
                  ncol  = length(unique(res_table$Prog)),
                  byrow = T)
  
  rownames(fdr_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(fdr_mat)<-colnames(inDf)
  
  return(list(table = res_table,
              r2 = r2_mat,
              P = p_mat,
              FDR = fdr_mat))
}





my_adonis_terms_adjAbun<-function(distList, inDf, covDf, covar,info){
  #distList<-PRJEB22863_msv_dist_std
  #inDf<-covar_PRJEB22863
  #covDf<-PRJEB22863_covar
  #covar<-covar2
  #info<-info
  
  res_table<-NULL
  for (i in 1:length(distList)) { # 
    #i<-8
    cat(paste(i,"\n")) 
    
    n<-names(distList)[i] %>% 
      str_replace_all("msv_","") %>% 
      match(., info$organism)
    
    inCovar<-c(covar, info$organism[n])

    if(is.na(distList[[i]])==F){
      if (nrow(distList[[i]])>10 & ncol(distList[[i]])>10){
         adonis_res<-my_adonis_terms(distList[[i]], inDf, covDf = covDf, covar = inCovar)
          adonis_res_df<-data.frame(Species = rep(info$Short_name[n], dim(adonis_res)[1]),
                              Prog = colnames(inDf),
                              as.data.frame(adonis_res))
    
        res_table<-rbind(res_table, adonis_res_df)
      }
      else {
        adonis_res_df<-data.frame(Species = rep(info$Short_name[n], ncol(inDf)),
                              Prog = colnames(inDf),
                              R2=NA,P=NA,N=NA)
        res_table<-rbind(res_table, adonis_res_df)
       }
     }
    else {
        adonis_res_df<-data.frame(Species = rep(info$Short_name[n], ncol(inDf)),
                              Prog = colnames(inDf),
                              R2=NA,P=NA,N=NA)
        res_table<-rbind(res_table, adonis_res_df)
     }
  }
  
  res_table$Species<-as.character(res_table$Species)
  res_table$Prog<-as.character(res_table$Prog)
  
  res_table$fdr<-p.adjust(res_table$P,method = 'fdr')
  res_table$bonferroni<-p.adjust(res_table$P,method = 'bonferroni')
  
  # R2
  r2_mat<-matrix(data  = res_table$R2,
                 nrow  = length(unique(res_table$Species)),
                 ncol  = length(unique(res_table$Prog)),
                 byrow = T)
  
  rownames(r2_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(r2_mat)<-colnames(inDf)
  
  # P
  p_mat<-matrix(data  = res_table$P,
                nrow  = length(unique(res_table$Species)),
                ncol  = length(unique(res_table$Prog)),
                byrow = T)
  
  rownames(p_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(p_mat)<-colnames(inDf)
  
  # fdr
  fdr_mat<-matrix(data  = res_table$fdr,
                  nrow  = length(unique(res_table$Species)),
                  ncol  = length(unique(res_table$Prog)),
                  byrow = T)
  
  rownames(fdr_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(fdr_mat)<-colnames(inDf)
  
  return(list(table = res_table,
              r2 = r2_mat,
              P = p_mat,
              FDR = fdr_mat))
}



my_meta_p <- function(inVec, study_name, n_col, p_col) {
  require(metap)
  #require(metafor)
  
  #inVec<-cbind_adonis_res.table[20,]
  #study_name<-c("LLD", "300OB")
  #n_col<-c(10,15)
  #_col<-c(9,14)
  
  study_n <- inVec[n_col] %>% as.numeric
  study_p <- inVec[p_col] %>% as.numeric
  
  inDf <- data.frame(study_name, study_n, study_p)
  inDf<-subset(inDf,is.na(study_p)==F,drop=T)

  study_n<-na.omit(study_n)
  study_p<-na.omit(study_p)
  
  Wi=sqrt(study_n)
  #Convert p-values to Z-scores
  Zi=qnorm(1-(study_p/2))
  # Meta-zscore  
  Z=(sum(Zi*Wi)/sqrt(sum(Wi^2)))
  # Convert Z-score to p-value
  MetaP=2*pnorm(-abs(Z))
  
  #Cochran Q-test (test heterogeneity of your meta-analysis)
  
  WeightedZ= sum(sqrt(study_n)*Zi)
  totalSample=sum(study_n)
  
  #Calculate expected Z
  expected_Z= sqrt(study_n)*WeightedZ/totalSample
  het_Sum= sum((Zi-expected_Z)*(Zi-expected_Z))
  
  #Get p-value of heterogeneity test!
  my_pvalue_co=pchisq(het_Sum, lower.tail = F, df=length(study_p)-1)
  
  return(
    list(
      Meta.p = MetaP,
      Meta.hetero = het_Sum,
      Meta.hetero.p = my_pvalue_co
    )
  )
}

## Batch meta-analysis
my_batch_meta_p <- function(inDf,study_name,n_col,p_col,row_var_col = 1,col_var_col = 2) {
  #inDf<-cbind_adonis_res.table[c(1:10),]
  #study_name<-c("LLD", "300OB")
  #n_col<-c(10,15)
  #p_col<-c(9,14)
  #row_var_col<-1
  #col_var_col<-2
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_p,
    study_name = study_name,
    n_col = n_col,
    p_col = p_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 3, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.p',  "Meta.hetero", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  rowName <- unique(inDf[, row_var_col])
  colName <- unique(inDf[, col_var_col])
  
  N_row <- length(rowName)
  N_col <- length(colName)
  
  
  batch_res.p <- matrix(batch_res.unlist[, 1], ncol = N_col, byrow = F)
  colnames(batch_res.p) <- colName
  rownames(batch_res.p) <- rowName
  batch_res.p[is.na(batch_res.p)] <- 0
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  # fdr matrix
  batch_res.fdr <- matrix(
    data  = batch_res_edge$Meta.fdr.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = T
  )
  colnames(batch_res.fdr) <- colName
  rownames(batch_res.fdr) <- rowName
  batch_res.fdr[is.na(batch_res.fdr)] <- 1
  
  # bonferroni matrix
  batch_res.bon <- matrix(
    data  = batch_res_edge$Meta.bonferroni.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = T
  )
  
  colnames(batch_res.bon) <- colName
  rownames(batch_res.bon) <- rowName
  batch_res.bon[is.na(batch_res.bon)] <- 1
  
  return(
    list(
      table      = batch_res_edge,
      p          = batch_res.p,
      fdr        = batch_res.fdr,
      bonferroni = batch_res.bon
    )
  )
}



## Batch meta-analysis
lr_batch_meta_p <- function(inDf,study_name,n_col,p_col,row_var_col,prog) {
  #inDf<-cbind_dsv_res.table
  #study_name<-c("PRJEB22863","PRJNA397906","PRJNA751792","PRJNA762360","PRJNA770295","PRJEB43119")
  #n_col<-seq(3,13,by=2)
  #p_col<-seq(2,13,by=2)
  #row_var_col<-1
  #prog<-"response"
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_p,
    study_name = study_name,
    n_col = n_col,
    p_col = p_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 3, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.p',  "Meta.hetero", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  rowName <- unique(inDf[, row_var_col])
  
  N_row <- length(rowName)
  N_col <- 1
  
  batch_res.p <- matrix(batch_res.unlist[, 1], ncol = N_col, byrow = F)
  colnames(batch_res.p) <- prog
  rownames(batch_res.p) <- rowName
  batch_res.p[is.na(batch_res.p)] <- 0
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  # fdr matrix
  batch_res.fdr <- matrix(
    data  = batch_res_edge$Meta.fdr.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = T
  )
  colnames(batch_res.fdr) <- prog
  rownames(batch_res.fdr) <- rowName
  batch_res.fdr[is.na(batch_res.fdr)] <- 1
  
  # bonferroni matrix
  batch_res.bon <- matrix(
    data  = batch_res_edge$Meta.bonferroni.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = T
  )
  
  colnames(batch_res.bon) <- prog
  rownames(batch_res.bon) <- rowName
  batch_res.bon[is.na(batch_res.bon)] <- 1
  
  return(
    list(
      table      = batch_res_edge,
      p          = batch_res.p,
      fdr        = batch_res.fdr,
      bonferroni = batch_res.bon
    )
  )
}






## The Normal Quantile Transformation
qtrans<-function(x){
  k<-!is.na(x)
  k<-which(x!="-999")
  ran<-rank(as.numeric(x[k]))
  y<-qnorm((1:length(k)-0.5)/length(k))
  x[k]<-y[ran]
  x
}



## linear model
lm_btw_mats<-function(y_mat,x_mat,cov_mat,covar){
  require(reshape2)
  require(R.utils)
  
  ## test block   
  #  y_mat<- all_ba                 
  #  x_mat<- all_abun_clr[,-86]     
  #  cov_mat<- all_covar           
  #  covar<- covar1               
  ## test block
  
  my_lm<-function(y,x){
   # y<-y_mat[,1]
   # x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    #y_uniq_N        <- NA
    #x_uniq_N        <- NA
    #y_non_zero_N    <-NA
    #x_non_zero_N    <-NA
    #y_non_zero_rate <-NA
    #x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    #lm_input<-data.frame(Y = y, X = x) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    #y_uniq_N   <- length(unique(lm_input[,"Y"]))
    #x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    #y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    #x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    #y_non_zero_rate<-y_non_zero_N/N
    #x_non_zero_rate<-x_non_zero_N/N
    
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    try(lm_res <- summary(lm(Y~.,data = lm_input)), silent = T)
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N)),
                    #y_uniq_N     = y_uniq_N,
                    #x_uniq_N     = x_uniq_N,
                    #y_non_zero_N = y_non_zero_N,
                    #x_non_zero_N = x_non_zero_N,
                    #y_non_zero_rate = y_non_zero_rate,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 4, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}





## Meta-analysis
my_meta_lm <- function(inVec, study_name, beta_col, se_col) {
  require(meta)
  require(metafor)
  
  ## test start
  #inVec<-inDf[80,]
  #study_name<-c("PRJEB22863", "PRJNA397906","PRJNA399742","PRJNA751792","PRJNA762360","PRJNA770295")
  #beta_col<-seq(15,86,by=12)
  #se_col<-seq(16,86,by=12)
  ## test end
  
  study_beta <- inVec[beta_col] %>% as.numeric
  study_se <- inVec[se_col] %>% as.numeric
  
  #study_beta<-c(0.7, 0.65)
  #study_se<-c(0.07, 0.08)
  #stydy_n<-c(1000, 300)
  
  in_file <- data.frame(study_name, study_beta, study_se)
  
  m.hksj <- metagen(
    study_beta,
    study_se,
    data = in_file,
    studlab = study_name
  )
  
  m.hksj.res <- summary(m.hksj)
  
  return(
    list(
      Meta.beta = m.hksj.res$random$TE,
      Meta.se = m.hksj.res$random$seTE,
      Meta.p = m.hksj.res$random$p,
      Meta.I2 = m.hksj.res$I2,
      Meta.hetero.p = m.hksj$pval.Q
    )
  )
}




## Batch meta-analysis
my_batch_meta_lm <- function(inDf,study_name,beta_col,se_col) {
  ## test start
  #inDf<-vsv_crass_lm_res.edge[c(1:10),]
  #study_name<-c("LLD", "300OB","IBD")
  #beta_col<-c(3,10,17)
  #se_col<-c(4,11,18)
  ## test end
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_lm,
    study_name = study_name,
    beta_col = beta_col,
    se_col = se_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 5, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.beta', 'Meta.SE', "Meta.p", "Meta.I2", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  return(batch_res_edge)
}


lm_btw_mats_adjAbun<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col){
  #y_mat<-ob_contin[,c(1:10)]
  #x_mat<-ob_vsv[,c(1:100)]
  #cov_mat<-ob_basic
  #covar<-covar
  #abun<-ob_abun_clr
  #info<-info
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lm_btw_mats_adjAbun2<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col){
  #y_mat<-vsgv_lld
  #x_mat<-lld_exp[,c(1:10)]
  #cov_mat<-lld_basic
  #covar<-covar
  #abun<-lld_abun_clr
  #info<-info
  #abun_col<-9
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    i<-1
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]
    if(dim(y_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}


lm_btw_mats_adjAbunPCs<-function(y_mat, x_mat, cov_mat, covar, abun, pc_mat, info, abun_col){
  #  y_mat<-lld_ba[,c(1:10)]
  #  x_mat<-vsgv_lld[,c(1:100)]
  #  cov_mat<-lld_basic
  #  covar<-covar
  #  abun<-lld_abun_clr
  #  pc_mat<-lld_msv_pc_cum0.6
  #  info<-info
  
  cov_mat<-cbind(cov_mat, abun,pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
#    i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))] 
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lm_btw_mats_adjAbunPCs2<-function(y_mat, x_mat, cov_mat, covar, abun, pc_mat, info, abun_col){
#    y_mat<-vsgv_lld[,c(1:10)]
#    x_mat<-lld_exp[,c(1:10)]
#    cov_mat<-lld_basic
#    covar<-covar
#    abun<-lld_abun_clr
#    pc_mat<-lld_msv_pc_cum0.6
#    info<-info
#    abun_col<-9
  
  cov_mat<-cbind(cov_mat, abun,pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
#    i<-2
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]%>% as.matrix
    colnames(y_mat_i)<-colnames(y_mat)[grep(spe_name,colnames(y_mat))]
    
    if(dim(y_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

## Get density of 2-demision dataset
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


## linear model
lr_btw_mats<-function(y_mat,x_mat,cov_mat,covar, direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  # direction: 1 means samples in row and variables in column; 2 means samples in column and variables in row
  
  require(reshape2)
  require(glmnet)
  require(R.utils)
  
  ## test block
  #y_mat <- lld_disea[,c(1:5)]
  #x_mat <- lld_vsv[,c(1:3)]
  #cov_mat <- lld_covar
  #covar   <- covar #covar
  #direction<-c(1,1)
  ## test block
  
  my_lm<-function(y,x){
    #y<-y_mat[,2]
    #x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    y_0_N<-NA
    y_1_N<-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Gender"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"Gender"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    y_0_N<-sum(lm_input[,"Y"]==0)
    y_1_N<-sum(lm_input[,"Y"]==1)
    
    lm_input_tmp <- apply(lm_input[,c(2:ncol(lm_input))], 2, qtrans) %>% as.data.frame
    lm_input<-data.frame(Y = lm_input[,1],lm_input_tmp)
    
    try(lm_res <- summary(glm(Y~., data = lm_input, family = 'binomial')))
    
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate,
                    y_0_N = y_0_N,
                    y_1_N = y_1_N)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 12, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2:1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","y_0_N","y_1_N","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}


lr_btw_mats_adjAbun<-function(y_mat, x_mat, cov_mat, covar, abun,info, abun_col){
  #y_mat<-lld_intri[,c(1:10)]
  #x_mat<-lld_vsv[,c(1:100)]
  #cov_mat<-lld_basic
  #covar<-covar
  #abun<-lld_abun
  #info<-info
  #direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lr_btw_mats_adjAbun2<-function(y_mat, x_mat, cov_mat, covar, abun,info, abun_col){
  #y_mat<-lld_intri[,c(1:10)]
  #x_mat<-lld_vsv[,c(1:100)]
  #cov_mat<-lld_basic
  #covar<-covar
  #abun<-lld_abun
  #info<-info
  #direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]
    if(dim(y_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}


lr_btw_mats_adjAbunPCs<-function(y_mat, x_mat, cov_mat, covar, abun,pc_mat,info, abun_col){
 # y_mat<-lld_intri[,c(1:10)]
#  x_mat<-lld_vsv[,c(1:100)]
# cov_mat<-lld_basic
#  covar<-covar
#  abun<-lld_abun_clr
#  pc_mat<-
#  info<-info
#  direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun, pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lr_btw_mats_adjAbunPCs2<-function(y_mat, x_mat, cov_mat, covar, abun,pc_mat,info, abun_col){
  # y_mat<-lld_intri[,c(1:10)]
  #  x_mat<-lld_vsv[,c(1:100)]
  # cov_mat<-lld_basic
  #  covar<-covar
  #  abun<-lld_abun_clr
  #  pc_mat<-
  #  info<-info
  #  direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun, pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]%>% as.matrix
    colnames(y_mat_i)<-colnames(y_mat)[grep(spe_name,colnames(y_mat))]
    
    if(dim(y_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}



## Mediation analysis for linear model
my_lm_mediation<-function(input.inv,input.med, input.dv,  covDf){
  #input.inv<-lld_exp$melatonine
  #input.med<-lld_vsv$`Faecalibacterium cf. prausnitzii KLE1255:1373_1377`
  #input.dv <-lld_ba$CA_dehydro_deconju_ratio
  #covDf<-lld_covar[,c('Gender','Age','BMI','Reads_number','s__Faecalibacterium_prausnitzii')]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=lm(input.dv~.,input.df.rmna[,-3])
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=lm(input.med~., input.df.rmna[,-2]) # input.inv
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=lm(input.dv~.,input.df.rmna)
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
    }
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
  }
  
  return(res_list)
}

## Bidirectional mediation analysis for linear model
lm_bimediation<-function(inVec, indvDf, dvDf1, dvDf2, covDf, covar ){
  ## test block
  #inVec<-exp_vsv_ba_df[1,]
  #indvDf<-lld_exp
  #dvDf1<-lld_vsv
  #dvDf2<-lld_ba
  #covDf<-lld_covar
  #covar<-c('Gender','Age','BMI','Reads_number')
  # test block
  
  if(is.na(inVec[4])){
    covar<-covar
  }else{
    covar<-c(covar, colnames(covDf)[grep(inVec[4],colnames(covDf))])
  }
  
  indv<-indvDf[,match(inVec[1], colnames(indvDf))]
  dv1<-dvDf1[,match(inVec[2], colnames(dvDf1))]
  dv2<-dvDf2[,match(inVec[3], colnames(dvDf2))]
  
  dir1_res <- my_lm_mediation(indv, dv1, dv2, covDf[,covar])
  dir2_res <- my_lm_mediation(indv, dv2, dv1, covDf[,covar])
  
  names(dir1_res)<-paste("dir1.",names(dir1_res),sep = "")
  names(dir2_res)<-paste("dir2.",names(dir2_res),sep = "")
  
  MediationDirection<-"none"
  if(!is.na(dir1_res$dir1.Prop.mediated.p) & !is.na(dir2_res$dir2.Prop.mediated.p)){
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "both"}
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p>0.05){MediationDirection <- "indv_dv1_dv2"}
    if( dir1_res$dir1.Prop.mediated.p>0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "indv_dv2_dv1"}
  }
  
  bires<-list(indv=inVec[1], dv1=inVec[2], dv2=inVec[3],MediationDirection = MediationDirection)
  
  res<-c(bires,dir1_res,dir2_res)
  
  return(res)
}


## Mediation analysis for logistic model
my_lr_mediation_1<-function(inpu.inv, input.dv, input.med, covDf){
  #input.inv<-indv
  #input.med<-dv2
  #input.dv <-dv1
  #covDf<- covDf[,covar]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  med<-input.df.rmna$input.med
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  input.df.rmna$input.med<-med
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=lm(input.dv~.,input.df.rmna[,-3])
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=glm(input.med~., input.df.rmna[,-2], family = "binomial")
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=lm(input.dv~.,input.df.rmna) 
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
      
    }
    
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
    
  }
  
  return(res_list)
}


my_lr_mediation_2<-function(inpu.inv, input.dv, input.med, covDf){
  #input.inv<-lld_exp$cereals
  #input.med<-lld_ba$GDCA
  #input.dv <-lld_dsv$`[Eubacterium] hallii DSM 3353:2969_2983`
  #covDf<-lld_covar[,c('Gender','Age','BMI','Reads_number','s__Eubacterium_hallii')]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  dv<-input.df.rmna$input.dv
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  input.df.rmna$input.dv<-dv
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=glm(input.dv~.,input.df.rmna[,-3], family = "binomial")
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=lm(input.med~., input.df.rmna[,-2])
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=glm(input.dv~.,input.df.rmna, family = "binomial") 
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
      
    }
    
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
    
  }
  
  return(res_list)
}





lr_bimediation<-function(inVec, indvDf, dvDf1, dvDf2, covDf, covar ){
  ## test block
  #inVec<-exp_dsv_ba_df[1,]
  #indvDf<-lld_exp
  #dvDf1<-lld_dsv
  #dvDf2<-lld_ba
  #covDf<-lld_covar
  #covar<-c('Gender','Age','BMI','Reads_number')
  # test block
  
  if(is.na(inVec[4])){
    covar<-covar
  }else{
    covar<-c(covar, inVec[4])
  }
  
  indv<-indvDf[,match(inVec[1], colnames(indvDf))]
  dv1<-dvDf1[,match(inVec[2], colnames(dvDf1))]
  dv2<-dvDf2[,match(inVec[3], colnames(dvDf2))]
  
  dir1_res <- my_lr_mediation_1(indv, dv2, dv1, covDf[,covar])
  dir2_res <- my_lr_mediation_2(indv, dv1, dv2, covDf[,covar])
  
  names(dir1_res)<-paste("dir1.",names(dir1_res),sep = "")
  names(dir2_res)<-paste("dir2.",names(dir2_res),sep = "")
  
  MediationDirection<-"none"
  if(!is.na(dir1_res$dir1.Prop.mediated.p) & !is.na(dir2_res$dir2.Prop.mediated.p)){
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "both"}
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p>0.05){MediationDirection <- "indv_dv1_dv2"}
    if( dir1_res$dir1.Prop.mediated.p>0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "indv_dv2_dv1"}
  }
  
  bires<-list(indv=inVec[1], dv1=inVec[2], dv2=inVec[3],MediationDirection = MediationDirection)
  
  res<-c(bires,dir1_res,dir2_res)
  
  return(res)
}


get_PCs<-function(inDist, eig.cutoff = 0.6){
  #inDist <- all_msv_dist_std[[1]]
  #eig.cutoff <- 0.6
  #eig.cutoff<-as.numeric(eig.cutoff)
  
  pcoa_res <- cmdscale(inDist, k=20, eig = T)
  
  total_eig <- 0
  i <- 0
  
  while (total_eig<eig.cutoff & i<10) {
    i<-i+1
    total_eig <- total_eig + pcoa_res$eig[i]/100
  }
  
  if(i<=1){
    pcoa<-data.frame(X1 = pcoa_res$points[, 1])
  }else{
    pcoa <- data.frame(pcoa_res$points)[,c(1:i)]
  }
  
  return(list(PCoA= pcoa,PC_num = i,Total_eig = total_eig))
  
}


my_uniq<-function(x){
  length(unique(x))
}



## Recalculate gene position
reCalcuPos<-function(spe, spe_scaf){
  spe_gene<-spe
  
  spe_gene<-data.frame(Taxonomy_id = str_replace_all(spe$contig_id, "\\..*", ""),
                       Contig_id = str_replace_all(spe$contig_id, ".*\\.", ""),
                       spe)
  Base_length<-NULL
  Base_length[1]<-0
  con_id<-spe_gene[1,2]
  for (i in 1:nrow(spe_gene)) {
    if(spe_gene[i,7] > spe_gene[i, 8]){
      tmp<-spe_gene[i,7]
      spe_gene[i,7]<-spe_gene[i, 8]
      spe_gene[i,8]<-tmp
    }
    
    if(i > 1){
      if( spe_gene[i,2]==con_id){
        Base_length[i]<-Base_length[i-1]
      }else{
        Base_length[i]<-spe_scaf$Length[match(con_id,spe_scaf$Contig)]+Base_length[i-1]
        con_id<-spe_gene[i,2]
      }
    }
  }
  
  spe_gene<-data.frame(spe_gene,
                       Base_length = Base_length,
                       Total_start = spe_gene$start + Base_length,
                       Total_stop  = spe_gene$stop  + Base_length)
  
  
}




## Clustering analysis
my_cluster<-function(inDist, ps.cutoff=0.8,my_seed = 1){
  require(NbClust)
  require(fpc)
  require(tsne)
  require(ggsci)
  
  # test
  #inDist<-all_msv_dist_std[[5]]
  # test

  if(nrow(inDist)>50){
  
    clu_n<-prediction.strength(as.dist(inDist), Gmin=2, Gmax=10, M=50,
                             clustermethod=claraCBI ,usepam=T,diss = T,
                             classification="centroid",cutoff=ps.cutoff,
                             distances=T,count=F)
    clu<-claraCBI(as.dist(inDist), clu_n$optimalk, usepam=T,diss=T)
    clu_df<-as.data.frame(clu$partition)
  
    rownames(clu_df)<-rownames(inDist)
  
    pcoa_res<-cmdscale(inDist, k=5, eig = T)
    pcoa <- data.frame(pcoa_res$points)
  
    set.seed(my_seed)
    tsne_res<-Rtsne::Rtsne(inDist, is_distance = TRUE,
                         perplexity = as.integer((nrow(inDist)-1)/3),
                         theta = 0, pca = T,
                         eta=as.integer(nrow(inDist)/12))
  
    tsne_df <- tsne_res$Y %>%
      data.frame() %>%
      setNames(c("X", "Y"))
     tsne_df <- data.frame(Cluster = as.factor(clu$partition),tsne_df)
  
     return(list(pcoa=pcoa, pcoa_res = pcoa_res, clu_n=clu_n,tsne_df=tsne_df))
   }
   else{
     return(list(pcoa=NA, pcoa_res = NA, clu_n=NA,tsne_df=NA))
  }
}




## Permutation kruskal-wallis test
permKW_btw_mats<-function(mat0,mat1,direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  require(reshape2)
  require(rstatix)
  require(coin)
  
  #mat0<-all_prog
  #mat1<-all_msv_cluster_sub
  #direction<-c(1,1)
  
  if(direction[1]==2){mat0 =t(mat0)}
  if(direction[2]==2){mat1 =t(mat1)}
  col_var<-mat0
  row_var<-mat1
  
  my_kw<-function(a,b){
    #a<-col_var[,1]
    #b<-row_var[,3]
    
    effSize  <- NA
    p.value  <- NA
    perm.fdr <- NA
    N        <- NA
    N_minor  <- NA
    uniq_N   <- NA
    
    kw_input <-data.frame(phen = a,
                          mbio = b)
    kw_input<-na.omit(kw_input)

    a<-table(kw_input$phen)
    b<-table(kw_input$mbio)

    if(length(a)==1 | length(b)==1){
      effSize   <- NA
      p.value   <- NA
      perm.fdr  <- NA
      N         <- NA
      }
    else if (nrow(kw_input)<10){
      effSize   <- NA
      p.value   <- NA
      perm.fdr  <- NA
      N         <- NA
    }
    else if (nrow(kw_input)>=10){
      res.kw <- kruskal_test(phen ~ as.factor(mbio), data = kw_input)
      res.kw.eff<-kruskal_effsize(phen ~ as.factor(mbio), data = kw_input)
      res.kw.perm <- oneway_test(phen ~ as.factor(mbio), data=kw_input)
    
      effSize   <- res.kw.eff$effsize
      p.value   <- pvalue(res.kw)
      perm.fdr  <- pvalue(res.kw.perm)
      N         <- nrow(kw_input)
      }

    try(return(list(effSize = effSize,
                    p.value = p.value,
                    perm.fdr=perm.fdr, 
                    N=N)),silent = T)
  }
  
  col_var_row_var<-sapply(
    as.data.frame(col_var),
    function(x) Map(function(a,b) my_kw(a,b),
                    list(x),
                    as.data.frame(row_var)
    )
  )

  col_var_row_var.unlist <- matrix(unlist(col_var_row_var), ncol =4, byrow = T)
  
  # beta matrix
  col_var_row_var.beta<-matrix(col_var_row_var.unlist[,1],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.beta)<-colnames(col_var)
  rownames(col_var_row_var.beta)<-colnames(row_var)
  col_var_row_var.beta[is.na(col_var_row_var.beta)]<-0
  
  # p matrix
  col_var_row_var.p <- matrix(col_var_row_var.unlist[,2],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.p)<-colnames(col_var)
  rownames(col_var_row_var.p)<-colnames(row_var)
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  # convert matrix to edge list
  col_var_row_var_edge_p     <- melt(col_var_row_var.p)
  col_var_row_var_edge_beta  <- melt(col_var_row_var.beta)
  
  col_var_row_var_edge<-cbind(col_var_row_var_edge_beta,col_var_row_var.unlist)[,-c(3)]
  colnames(col_var_row_var_edge)<-c("Taxa", "Phenotype", "effectSize", "p","perm.fdr","N")
  
  # add p adjust
  #col_var_row_var_edge<-data.frame(as.data.frame(col_var_row_var_edge),
  #                                 fdr.p = p.adjust(col_var_row_var_edge$p, method = "fdr"),
  #                                 bonferroni.p = p.adjust(col_var_row_var_edge$p, method = "bonferroni"))
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var_edge$perm.fdr,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(col_var)
  rownames(col_var_row_var.fdr)<-colnames(row_var)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  
  return(list(table      = col_var_row_var_edge,
              effectSize = col_var_row_var.beta,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr))
  
}



logistic_btw_mats_noadjAbun_dsv<-function(y_mat, x_mat, cov_mat, covar, info, pheno){

  ############# test block  
  #y_mat<-mel_prog
  #x_mat<-mel_vsgv
  #cov_mat<-mel_basic
  #covar<-covar
  #info<-info
  #pheno<-"response_code"
  ############# test block

  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
        covar_i<-covar
        y_x_i <- logistic_btw_mats_dsv(y_mat, x_mat_i, cov_mat, covar_i,pheno)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  
  y_x.edge<-as.data.frame(y_x.edge)
  return(y_x.edge)
}



logistic_btw_mats_adjAbun_dsv<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col,pheno){

  ############# test block  
  #y_mat<-PRJNA762360_prog
  #x_mat<-PRJNA762360_dsgv
  #cov_mat<-PRJNA762360_basic
  #covar<-PRJNA762360_covar
  #abun<-PRJNA762360_abun_clr
  #info<-info
  #abun_col<-1
  #pheno<-"response_code"
  ############# test block

  cov_mat<-cbind(cov_mat, abun)
  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-55
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- logistic_btw_mats_dsv(y_mat, x_mat_i, cov_mat, covar_i,pheno)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- logistic_btw_mats_dsv(y_mat, x_mat_i, cov_mat, covar_i,pheno)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  return(y_x.edge)
}


## logistic regression model

logistic_btw_mats_dsv<-function(y_mat,x_mat,cov_mat,covar,pheno){
  require(reshape2)
  require(R.utils)
  
  ## test block   
  #  y_mat<- y_mat               
  #  x_mat<- x_mat_i    
  #  cov_mat<- cov_mat         
  #  covar<- covar_i
  #  pheno<-pheno       
  ## test block

  my_lg<-function(pheno,x){
   # y<-y_mat[,1]
   # x<-x_mat_i[,2]
     y<-y_mat[,pheno]
    
    or    <- NA
    se      <- NA
    p.value <- NA
    or_left<-NA
    or_right<-NA

    N       <- NA
    
    #x_uniq_N        <- NA
    #x_non_zero_N    <-NA
    #x_non_zero_rate <-NA
    
    lg_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lg_input)
    
    #x_uniq_N   <- length(unique(lg_input[,"X"]))
    #x_non_zero_N <- sum(!isZero(lg_input[,"X"]))
    #x_non_zero_rate<-x_non_zero_N/N
    
    ##lg_input <- apply(lg_input, 2, qtrans) %>% as.data.frame
    lg_input<-as.data.frame(lg_input)
    
    if(nrow(lg_input)< 10 | min(table(lg_input$Y))<3 | length(table(lg_input$Y))==1 | length(table(lg_input$X))==1){
        try(return(list(or =NA,
                    se =NA,
                    p.value =NA,
                    or_left=NA,
                    or_right=NA,
                    N = N)),
                    #x_uniq_N     = x_uniq_N
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
      }

   else{
    
    try(lg_res_model <- glm(Y~.,family=binomial(logit), data=lg_input), silent = T)
    lg_res <- summary(lg_res_model)

    indv<-'X'
    
    try(or    <- exp(lg_res$coefficients[match(indv,rownames(lg_res$coefficients)),1]), silent = T)
    try(se      <- lg_res$coefficients[match(indv,rownames(lg_res$coefficients)),2], silent = T)
    try(p.value <- lg_res$coefficients[match(indv,rownames(lg_res$coefficients)),4],silent = T)

    if(is.na(p.value)==F & (p.value>0.97 | p.value<2e-16)){
      or_left<-"NA"
      or_right<-"NA"
    }
    else{
      ci<-confint(lg_res_model)
      ci<-unlist(ci)
      try(or_left <- exp(ci[match(indv,rownames(ci)),1]),silent = T)
      try(or_right <- exp(ci[match(indv,rownames(ci)),2]),silent = T)
    }

    try(return(list(or = or,
                    se = se,
                    p.value = p.value,
                    or_left=or_left,
                    or_right=or_right,
                    N = N)),
                    #x_uniq_N     = x_uniq_N,
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    }
  }
  
  y_x<-matrix(,nrow=ncol(x_mat),ncol=6)

  for(i in 1:ncol(x_mat)){
    # i<-5
    #result<-unlist(my_cox("os","os_event",x_mat[,i]))
    result<-unlist(my_lg(pheno,x_mat[,i]))
    y_x[i,]<-matrix(result,ncol=6)
    }

  rownames(y_x)<-colnames(x_mat)
  
  y_x_edge<-data.frame(colnames(x_mat),as.data.frame(y_x),
                       fdr.p = p.adjust(y_x[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x[,3], method = "bonferroni"))
  
  colnames(y_x_edge)<-c("X", "OR","SE", "p","OR_left","OR_right","N","fdr.p","bonferroni.p")
  return(y_x_edge)

}



logistic_btw_mats_noadjAbun_vsv<-function(y_mat, x_mat, cov_mat, covar, info, pheno){

  ############# test block  
  #y_mat<-mel_prog
  #x_mat<-mel_vsgv
  #cov_mat<-mel_basic
  #covar<-covar
  #info<-info
  #pheno<-"response_code"
  ############# test block

  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
        covar_i<-covar
        y_x_i <- logistic_btw_mats_vsv(y_mat, x_mat_i, cov_mat, covar_i,pheno)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  
  y_x.edge<-as.data.frame(y_x.edge)
  return(y_x.edge)
}


logistic_btw_mats_adjAbun_vsv<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col,pheno){

  ############# test block  
  #y_mat<-mel_prog
  #x_mat<-mel_vsgv
  #cov_mat<-mel_basic
  #covar<-covar
  #abun<-mel_abun_clr 
  #info<-info
  #abun_col<-9
  #pheno<-"response_code"
  ############# test block

  cov_mat<-cbind(cov_mat, abun)
  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- logistic_btw_mats_vsv(y_mat, x_mat_i, cov_mat, covar_i,pheno)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- logistic_btw_mats_vsv(y_mat, x_mat_i, cov_mat, covar_i,pheno)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  return(y_x.edge)
}


## logistic regression model

logistic_btw_mats_vsv<-function(y_mat,x_mat,cov_mat,covar,pheno){
  require(reshape2)
  require(R.utils)
  
  ## test block   
  #  y_mat<- all_prog                 
  #  x_mat<- all_abun_clr[,-ncol(all_abun_clr)]     
  #  cov_mat<- all_covar           
  #  covar<- covar              
  ## test block
  
  my_lg<-function(pheno,x){
   # y<-y_mat[,1]
   # x<-x_mat[,i]
     y<-y_mat[,pheno]
    
    or    <- NA
    se      <- NA
    p.value <- NA
    or_left<-NA
    or_right<-NA

    N       <- NA
    
    #x_uniq_N        <- NA
    #x_non_zero_N    <-NA
    #x_non_zero_rate <-NA

    qt<-quantile(x, probs = c(0.25,0.5,0.75),na.rm=T)
  
     subtype<-0
     for (j in 1:length(x)){
      if (is.na(x[j])==T) {subtype[j]=NA}
      else if (x[j]<=qt[1]){subtype[j]=1}
      else if (x[j]>qt[1] & x[j]<=qt[2]){subtype[j]=2}
      else if (x[j]>qt[2] & x[j]<=qt[3]){subtype[j]=3}
      else if (x[j]>qt[3]){subtype[j]=4}
     }
    
    lg_input<-data.frame(Y = y, X = subtype,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lg_input)
    
    #x_uniq_N   <- length(unique(lg_input[,"X"]))
    #x_non_zero_N <- sum(!isZero(lg_input[,"X"]))
    #x_non_zero_rate<-x_non_zero_N/N
    
    ##lg_input <- apply(lg_input, 2, qtrans) %>% as.data.frame
    lg_input<-as.data.frame(lg_input)
    
    if(nrow(lg_input)<10 | min(table(lg_input$Y))<3){
        try(return(list(or =NA,
                    se =NA,
                    p.value =NA,
                    or_left=NA,
                    or_right=NA,
                    N = N)),
                    #x_uniq_N     = x_uniq_N,
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
      }
    
    try(lg_res_model <- glm(Y~.,family=binomial(logit), data=lg_input), silent = T)
    lg_res <- summary(lg_res_model)

    indv<-'X'
    
    try(or    <- exp(lg_res$coefficients[match(indv,rownames(lg_res$coefficients)),1]), silent = T)
    try(se      <- lg_res$coefficients[match(indv,rownames(lg_res$coefficients)),2], silent = T)
    try(p.value <- lg_res$coefficients[match(indv,rownames(lg_res$coefficients)),4],silent = T)

    if(is.na(p.value)==F & (p.value>0.97 | p.value<2e-16)){
      or_left<-"NA"
      or_right<-"NA"
    }
    else{
      ci<-confint(lg_res_model)
      ci<-unlist(ci)
      try(or_left <- exp(ci[match(indv,rownames(ci)),1]),silent = T)
      try(or_right <- exp(ci[match(indv,rownames(ci)),2]),silent = T)
    }

    try(return(list(or = or,
                    se = se,
                    p.value = p.value,
                    or_left=or_left,
                    or_right=or_right,
                    N = N)),
                    #x_uniq_N     = x_uniq_N,
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  y_x<-matrix(,nrow=ncol(x_mat),ncol=6)

  for(i in 1:ncol(x_mat)){
    # i<-5
    #result<-unlist(my_cox("os","os_event",x_mat[,i]))
    result<-unlist(my_lg(pheno,x_mat[,i]))
    y_x[i,]<-matrix(result,ncol=6)
    }

  rownames(y_x)<-colnames(x_mat)
  
  y_x_edge<-data.frame(colnames(x_mat),as.data.frame(y_x),
                       fdr.p = p.adjust(y_x[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x[,3], method = "bonferroni"))
  
  colnames(y_x_edge)<-c("X", "OR","SE", "p","OR_left","OR_right","N","fdr.p","bonferroni.p")
  return(y_x_edge)

}


logistic_btw_mats_adjAbunPCs<-function(y_mat, x_mat, cov_mat, covar, abun, pc_mat,info, abun_col,pheno){

  ############# test block  
  #y_mat<-all_prog
  #x_mat<-vsgv
  #cov_mat<-all_basic
  #covar<-covar
  #abun<-all_abun_clr
  #pc_mat<-msv_pc_cum0.6
  #info<-info
  #abun_col<-1
  #pheno<-"response_code"
  ############# test block

  cov_mat<-cbind(cov_mat,abun,pc_mat)
  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]

    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- logistic_btw_mats(y_mat, x_mat_i, cov_mat, covar_i,pheno)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- logistic_btw_mats(y_mat, x_mat_i, cov_mat, covar_i,pheno)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  return(y_x.edge)
}



cox_btw_mats_adjAbunPCs<-function(y_mat, x_mat, cov_mat, covar, abun, pc_mat,info, abun_col,surv_time,surv_event){

  ############# test block  
  #  y_mat<-all_surv
  #  x_mat<-vsgv[,c(1:100)]
  #  cov_mat<-all_covar
  #  covar<-covar
  #  abun<-all_abun_clr[,-ncol(all_abun_clr)] 
  #  info<-info
  #  abun_col<-9
  ############# test block

  cov_mat<-cbind(cov_mat,abun,pc_mat)
  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- cox_btw_mats(y_mat, x_mat_i, cov_mat, covar_i,surv_time,surv_event)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- cox_btw_mats(y_mat, x_mat_i, cov_mat, covar_i,surv_time,surv_event)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  #y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  #y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  return(y_x.edge)
}




cox_btw_mats_noadjAbun_dsv<-function(y_mat, x_mat, cov_mat, covar, info,surv_time,surv_event){

  ############# test block  
  #  y_mat<-mel_surv
  #  x_mat<-mel_vsgv
  #  cov_mat<-mel_basic
  #  covar<-mel_covar
  #  info<-info
  #  surv_time<-"os"
  #  surv_event<-"os_event"
  ############# test block

  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
        covar_i<-covar
        y_x_i <- cox_btw_mats_dsv(y_mat, x_mat_i, cov_mat, covar_i,surv_time,surv_event)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  #y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  #y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  return(y_x.edge)
}




cox_btw_mats_adjAbun_dsv<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col,surv_time,surv_event){

  ############# test block  
  #  y_mat<-PRJNA762360_surv
  #  x_mat<-PRJNA762360_dsgv
  #  cov_mat<-PRJNA762360_basic
  #  covar<-PRJNA762360_covar
  #  abun<-PRJNA762360_abun_clr
  #  info<-USA_UK_info
  #  abun_col<-1
  #  surv_time<-"os"
  #  surv_event<-"os_event"
  ############# test block

  cov_mat<-cbind(cov_mat, abun)
  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- cox_btw_mats_dsv(y_mat, x_mat_i, cov_mat, covar_i,surv_time,surv_event)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- cox_btw_mats_dsv(y_mat, x_mat_i, cov_mat, covar_i,surv_time,surv_event)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  #y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  #y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  return(y_x.edge)
}


cox_btw_mats_dsv<-function(y_mat,x_mat,cov_mat,covar,surv_time,surv_event){
  require(reshape2)
  require(R.utils)
  
  ## test block   
  #  y_mat<- y_mat                
  #  x_mat<- x_mat_i    
  #  cov_mat<- cov_mat          
  #  covar<- covar_i
  #  surv_time<-"os"
  #  surv_event<-"os_event"
  ## test block

  
  my_cox<-function(time,event,x){

   surv_time<-y_mat[,time]
   surv_event<-y_mat[,event]
   
   # surv_time<-y_mat[,"os"]
   # surv_event<-y_mat[,"os_event"]
   # x<-x_mat[,1]
    
    hr    <- NA
    se      <- NA
    p.value <- NA
    hr_left<-NA
    hr_right<-NA
    N       <- NA
    #x_uniq_N        <- NA
    #x_non_zero_N    <-NA
    #x_non_zero_rate <-NA
    
    cox_input<-data.frame(surv_time =surv_time,surv_event=surv_event,X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(cox_input)
    
    #x_uniq_N   <- length(unique(cox_input[,"X"]))
    #x_non_zero_N <- sum(!isZero(cox_input[,"X"]))
    #x_non_zero_rate<-x_non_zero_N/N
    
    ##cox_input <- apply(lg_input, 2, qtrans) %>% as.data.frame
    
    cox_input<-as.data.frame(cox_input)
     cox_input<-subset(cox_input,is.na(X)!=T & is.na(surv_time)!=T & is.na(surv_event)!=T, drop=T)
     if(is.null(nrow(cox_input))==T){
        try(return(list(hr =NA,
                    se = NA,
                    p.value = NA,
                    hr_left=NA,
                    hr_right=NA,
                    N = N)),
                    #x_uniq_N     = x_uniq_N,
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
      }
      else if (is.null(nrow(cox_input))==F & nrow(cox_input)<20){
        try(return(list(hr =NA,
                    se = NA,
                    p.value = NA,
                    hr_left=NA,
                    hr_right=NA,
                    N = N)),
                    #x_uniq_N     = x_uniq_N)),
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
      }
      else{
       try(cox_res_model <- coxph(Surv(surv_time, surv_event)~.,cox_input,na.action=na.exclude), silent = T)
       cox_res <- summary(cox_res_model)

      indv<-'X'
    
      try(hr    <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),2], silent = T)
      try(se      <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),3], silent = T)
      try(p.value <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),5],silent = T)
      try(hr_left <- cox_res$conf.int[match(indv,rownames(cox_res$conf.int)),3],silent = T)
      try(hr_right <- cox_res$conf.int[match(indv,rownames(cox_res$conf.int)),4],silent = T)

      try(return(list(hr = hr,
                    se = se,
                    p.value = p.value,
                    hr_left=hr_left,
                    hr_right=hr_right,
                    N = N)),
                    #x_uniq_N     = x_uniq_N
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    }
  }

  y_x<-matrix(,nrow=ncol(x_mat),ncol=6)

  for(j in 1:ncol(x_mat)){
    # i<-5
    #result<-unlist(my_cox("os","os_event",x_mat[,j]))
    result<-unlist(my_cox(surv_time,surv_event,x_mat[,j]))
    y_x[j,]<-matrix(result,ncol=6)
    }

  rownames(y_x)<-colnames(x_mat)
  
  y_x_edge<-data.frame(colnames(x_mat),as.data.frame(y_x),
                       fdr.p = p.adjust(y_x[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x[,3], method = "bonferroni"))
  
  colnames(y_x_edge)<-c("X", "HR","SE", "p","HR_left","HR_right","N","fdr.p","bonferroni.p")
  return(y_x_edge)
}



cox_btw_mats_noadjAbun_vsv<-function(y_mat, x_mat, cov_mat, covar, info,surv_time,surv_event){

  ############# test block  
  #  y_mat<-mel_surv
  #  x_mat<-mel_vsgv
  #  cov_mat<-mel_basic
  #  covar<-mel_covar
  #  info<-info
  #  surv_time<-"os"
  #  surv_event<-"os_event"
  ############# test block

  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
        covar_i<-covar
        y_x_i <- cox_btw_mats_vsv(y_mat, x_mat_i, cov_mat, covar_i,surv_time,surv_event)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  #y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  #y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  return(y_x.edge)
}


cox_btw_mats_adjAbun_vsv<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col,surv_time,surv_event){

  ############# test block  
  #  y_mat<-mel_surv
  #  x_mat<-mel_vsgv
  #  cov_mat<-mel_basic
  #  covar<-mel_covar
  #  abun<-mel_abun_clr 
  #  info<-info
  #  abun_col<-9
  #  surv_time<-"os"
  #  surv_event<-"os_event"
  ############# test block

  cov_mat<-cbind(cov_mat, abun)
  y_x.edge<-NULL

  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- cox_btw_mats_vsv(y_mat, x_mat_i, cov_mat, covar_i,surv_time,surv_event)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- cox_btw_mats_vsv(y_mat, x_mat_i, cov_mat, covar_i,surv_time,surv_event)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  #y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  #y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  return(y_x.edge)
}


## cox regression model

cox_btw_mats_vsv<-function(y_mat,x_mat,cov_mat,covar,surv_time,surv_event){
  require(reshape2)
  require(R.utils)
  
  ## test block   
  #  y_mat<- y_mat                
  #  x_mat<- x_mat_i    
  #  cov_mat<- cov_mat          
  #  covar<- covar
  #  surv_time<-"os"
  #  surv_event<-"os_event"
  ## test block

  
  my_cox<-function(time,event,x){

   surv_time<-y_mat[,time]
   surv_event<-y_mat[,event]
   
   # surv_time<-y_mat[,"os"]
   # surv_event<-y_mat[,"os_event"]
   # x<-x_mat[,1]
    
    hr    <- NA
    se      <- NA
    p.value <- NA
    hr_left<-NA
    hr_right<-NA

    N       <- NA
    
    #x_uniq_N        <- NA
    #x_non_zero_N    <-NA
    #x_non_zero_rate <-NA
    
     qt<-quantile(x, probs = c(0.25,0.5,0.75),na.rm=T)
  
     subtype<-0
     for (j in 1:length(x)){
      if (is.na(x[j])==T) {subtype[j]=NA}
      else if (x[j]<=qt[1]){subtype[j]=1}
      else if (x[j]>qt[1] & x[j]<=qt[2]){subtype[j]=2}
      else if (x[j]>qt[2] & x[j]<=qt[3]){subtype[j]=3}
      else if (x[j]>qt[3]){subtype[j]=4}
     }
    
    cox_input<-data.frame(surv_time =surv_time,surv_event=surv_event,X = subtype,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(cox_input)
    
    #x_uniq_N   <- length(unique(cox_input[,"X"]))
    #x_non_zero_N <- sum(!isZero(cox_input[,"X"]))
    #x_non_zero_rate<-x_non_zero_N/N
    
    ##cox_input <- apply(lg_input, 2, qtrans) %>% as.data.frame
    
    cox_input<-as.data.frame(cox_input)
     cox_input<-subset(cox_input,is.na(X)!=T & is.na(surv_time)!=T & is.na(surv_event)!=T, drop=T)
     if(is.null(nrow(cox_input))==T){
        try(return(list(hr =NA,
                    se = NA,
                    p.value = NA,
                    hr_left=NA,
                    hr_right=NA,
                    N = N)),
                    #x_uniq_N     = x_uniq_N)),
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
      }
      else if (is.null(nrow(cox_input))==F & nrow(cox_input)< 20 ){
        try(return(list(hr =NA,
                    se = NA,
                    p.value = NA,
                    hr_left=NA,
                    hr_right=NA,
                    N = N)),
                    #x_uniq_N     = x_uniq_N,
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
      }
      else{
       try(cox_res_model <- coxph(Surv(surv_time, surv_event)~.,cox_input,na.action=na.exclude), silent = T)
       cox_res <- summary(cox_res_model)

      indv<-'X'
    
      try(hr    <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),2], silent = T)
      try(se      <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),3], silent = T)
      try(p.value <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),5],silent = T)
      try(hr_left <- cox_res$conf.int[match(indv,rownames(cox_res$conf.int)),3],silent = T)
      try(hr_right <- cox_res$conf.int[match(indv,rownames(cox_res$conf.int)),4],silent = T)

      try(return(list(hr = hr,
                    se = se,
                    p.value = p.value,
                    hr_left=hr_left,
                    hr_right=hr_right,
                    N = N)),
                    #x_uniq_N     = x_uniq_N,
                    #x_non_zero_N = x_non_zero_N,
                    #x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    }
  }

  y_x<-matrix(,nrow=ncol(x_mat),ncol=6)

  for(i in 1:ncol(x_mat)){
    # i<-5
    #result<-unlist(my_cox("pfs","pfs_event",x_mat[,i]))
    result<-unlist(my_cox(surv_time,surv_event,x_mat[,i]))
    y_x[i,]<-matrix(result,ncol=6)
    }

  rownames(y_x)<-colnames(x_mat)
  
  y_x_edge<-data.frame(colnames(x_mat),as.data.frame(y_x),
                       fdr.p = p.adjust(y_x[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x[,3], method = "bonferroni"))
  
  colnames(y_x_edge)<-c("X", "HR","SE", "p","HR_left","HR_right","N","fdr.p","bonferroni.p")
  return(y_x_edge)
}




cluster_cox_os<-function(all_survival,cluster_file,outfile){

   #cluster_file<-NSCLC_msv_cluster_sub
   #all_survival<-NSCLC_survival

    hr    <- c()
    se      <- c()
    p.value <- c()
    hr_left<-c()
    hr_right<-c()
    N       <- c()
    
   for(i in 1:(ncol(cluster_file)-1)){
    
    clust<-cluster_file[,c(i,ncol(cluster_file))]
    colnames(clust)<-c("X","id")
    clust<-na.omit(clust)

    num<- nrow(clust)
    N[i]<-num
     if(num<10){
        hr[i]<-NA
        se[i]<-NA
        p.value[i]<-NA
        hr_left[i]<-NA
        hr_right[i]<-NA
      }
      else{
       cox_input<-data.frame(merge(all_survival,clust,by.x="id",by.y="id"))
       cox_input<-subset(cox_input,is.na(X)==F & is.na(os)==F & is.na(os_event)==F,drop=T)
       try(cox_res_model <- coxph(Surv(os,os_event)~X,cox_input,na.action=na.exclude), silent = T)
       cox_res <- summary(cox_res_model)

      indv<-'X'
      try(hr[i]    <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),2], silent = T)
      try(se[i]     <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),3], silent = T)
      #try(p.value[i] <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),5],silent = T)
      try(p.value[i] <- unlist(cox_res)$sctest.pvalue,silent = T)
      try(hr_left[i] <- cox_res$conf.int[match(indv,rownames(cox_res$conf.int)),3],silent = T)
      try(hr_right[i] <- cox_res$conf.int[match(indv,rownames(cox_res$conf.int)),4],silent = T)
    }
  }

  re<-cbind(hr,se,p.value,hr_left,hr_right,N)
  colnames(re)<-c("HR","SE", "p","HR_left","HR_right","N")
  rownames(re)<-colnames(cluster_file)[-ncol(cluster_file)]

  write.csv(re,outfile)
}



cluster_cox_pfs<-function(all_survival,cluster_file,outfile){

   #cluster_file<-NSCLC_msv_cluster_sub
   #all_survival<-NSCLC_survival

    hr    <- c()
    se      <- c()
    p.value <- c()
    hr_left<-c()
    hr_right<-c()
    N       <- c()
    
   for(i in 1:(ncol(cluster_file)-1)){
    
    clust<-cluster_file[,c(i,ncol(cluster_file))]
    colnames(clust)<-c("X","id")
    clust<-na.omit(clust)

    num<- nrow(clust)
    N[i]<-num
     if(num<10){
        hr[i]<-NA
        se[i]<-NA
        p.value[i]<-NA
        hr_left[i]<-NA
        hr_right[i]<-NA
      }
      else{
       cox_input<-data.frame(merge(all_survival,clust,by.x="id",by.y="id"))
       cox_input<-subset(cox_input,is.na(X)==F & is.na(pfs)==F & is.na(pfs_event)==F,drop=T)
       try(cox_res_model <- coxph(Surv(pfs,pfs_event)~X,cox_input,na.action=na.exclude), silent = T)
       cox_res <- summary(cox_res_model)

      indv<-'X'
      try(hr[i]    <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),2], silent = T)
      try(se[i]     <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),3], silent = T)
      #try(p.value[i] <- cox_res$coefficients[match(indv,rownames(cox_res$coefficients)),5],silent = T)
      try(p.value[i] <- unlist(cox_res)$sctest.pvalue,silent = T)
      try(hr_left[i] <- cox_res$conf.int[match(indv,rownames(cox_res$conf.int)),3],silent = T)
      try(hr_right[i] <- cox_res$conf.int[match(indv,rownames(cox_res$conf.int)),4],silent = T)
    }
  }

  re<-cbind(hr,se,p.value,hr_left,hr_right,N)
  colnames(re)<-c("HR","SE", "p","HR_left","HR_right","N")
  rownames(re)<-colnames(cluster_file)[-ncol(cluster_file)]

  write.csv(re,outfile)
}





## Batch meta-analysis
my_batch_meta_or <- function(inDf,study_name,or_col,se_col) {
  ## test start
  #inDf<-vsv_crass_lm_res.edge[c(1:10),]
  #study_name<-c("LLD", "300OB","IBD")
  #beta_col<-c(3,10,17)
  #se_col<-c(4,11,18)
  ## test end
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_or,
    study_name = study_name,
    or_col = or_col,
    se_col = se_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 5, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.OR', 'Meta.SE', "Meta.p", "Meta.I2", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  return(batch_res_edge)
}



## Meta-analysis
my_meta_or <- function(inVec, study_name, or_col, se_col) {
  require(meta)
  require(metafor)
  
  ## test start
  #inVec<-inDf[80,]
  #study_name<-c("A", "B")
  #beta_col<-seq(15,86,by=12)
  #se_col<-seq(16,86,by=12)
  ## test end
  
  study_or <- inVec[or_col] %>% as.numeric
  study_se <- inVec[se_col] %>% as.numeric
  
  #study_or<-c(0.7, 0.65)
  #study_se<-c(0.07, 0.08)
  #stydy_n<-c(1000, 300)
  
  in_file <- data.frame(study_name, study_or, study_se)
  
  m.hksj <- metagen(
    log(study_or),
    study_se,
    data = in_file,
    sm="OR",
    studlab = study_name
  )
  
  m.hksj.res <- summary(m.hksj)
  
  return(
    list(
      Meta.OR = exp(m.hksj.res$random$TE),
      Meta.se = m.hksj.res$random$seTE,
      Meta.p = m.hksj.res$random$p,
      Meta.I2 = m.hksj.res$I2,
      Meta.hetero.p = m.hksj$pval.Q
    )
  )
}




my_batch_meta_p_lr<-function(datain,effect_col,se_col,dataout){
    #data<-mel_cbind_sv_abun_resp
     data<-datain
     nc<-ncol(data)
     nr<-nrow(data)
     
     #effect<-data[,seq(3,50,12)]
     #se<-data[,seq(4,50,12)]

     effect<-data[,effect_col]
     se<-data[,se_col]

     effect_meta<-c()
     pvalue_meta<-c()
     cilb<-c()
     ciub<-c()
     I2<-c()
     het_p<-c()

#i=944

  for(i in 1:nr){
    ### 删掉取值为0的对数，即se为无穷大的
    site_effect<-effect[i,]
    ind1<-which(site_effect==0|site_effect=="Inf"|is.na(site_effect)==T|site_effect=="NaN") 
    site_se<-se[i,]
    ind2<-which(site_se==0|site_se=="Inf"|is.na(site_se)==T) 
    ind<-unique(c(ind1,ind2))
    if(length(ind)==0){
       site_beta<-as.numeric(site_effect)
       site_se<-as.numeric(site_se)
       res<- rma(measure="GEN", yi=t(site_beta), sei=t(site_se),method="EE",control=list(maxiter=10000))
       effect_meta[i]<-res$b
       pvalue_meta[i]<-res$pval
       cilb[i]<-res$ci.lb
       ciub[i]<-res$ci.ub
       I2[i]<-res$I2
       het_p[i]<-res$QEp
     }
    else if (length(ind)>0 & length(ind)<=((nc-2)/6-2)){  
       site_beta<-log(as.numeric(site_effect[-ind]))
       site_se<-as.numeric(site_se[-ind])
       res<- rma(measure="GEN", yi=t(site_beta), sei=t(site_se),method="EE",control=list(maxiter=10000))
       effect_meta[i]<-res$b
       pvalue_meta[i]<-res$pval
       cilb[i]<-res$ci.lb
       ciub[i]<-res$ci.ub
       I2[i]<-res$I2
       het_p[i]<-res$QEp
    }
   else if (length(ind)==(nc-2)/6 | length(ind)==(nc-2)/6-1){ 
      effect_meta[i]<-NA
      pvalue_meta[i]<-NA
      cilb[i]<-NA
      ciub[i]<-NA
      I2[i]<-NA
      het_p[i]<-NA
   } 
 }

  re<-cbind(data,effect_meta,pvalue_meta,cilb,ciub,I2,het_p)
  re$fdr<-p.adjust(re$pvalue_meta,method="fdr")
 #write.csv(re, file = "08.Microbial_GWAS/test_mete_pvalue_cut15.csv")
  write.csv(re,dataout,row.names=F)
}






combine_meta_p_lr<-function(datain,dataout){
    #data<-nsclc_dsv_resp_res.table
     data<-datain
     nc<-ncol(data)
     nr<-nrow(data)
     effect<-data[,seq(2,nc,4)]
     se<-data[,seq(3,nc,4)]
     
     dataset_p<-data[,seq(4,nc,4)]

     effect_meta<-c()
     pvalue_meta<-c()
     se_meta<-c()
     cilb<-c()
     ciub<-c()
     I2<-c()
     het_p<-c()
     p_count<-c()
     #i=944

  for(i in 1:nr){
    ### 删掉取值为0的对数，即se为无穷大的
    site_effect<-effect[i,]
    ind1<-which(site_effect==0|site_effect=="Inf"|is.na(site_effect)==T) 
    site_se<-se[i,]
    ind2<-which(site_se==0|site_se=="Inf"|is.na(site_se)==T) 
    ind<-unique(c(ind1,ind2))
    
    site_p<-dataset_p[i,]

    if(length(ind)==0){
       site_beta<-log(as.numeric(site_effect))
       site_se<-as.numeric(site_se)
       res<- rma(measure="GEN", yi=t(site_beta), sei=t(site_se),method="EE",control=list(maxiter=10000))
       effect_meta[i]<-exp(res$b)
       se_meta[i]<-exp(res$se)
       pvalue_meta[i]<-res$pval
       cilb[i]<-exp(res$ci.lb)
       ciub[i]<-exp(res$ci.ub)
       I2[i]<-res$I2
       het_p[i]<-res$QEp
       
       p_count[i]<-length(which(site_p<=0.2)) 
     }
    else if (length(ind)>0 & length(ind)<=((nc-1)/4-2)){  
       site_beta<-log(as.numeric(site_effect[-ind]))
       site_se<-as.numeric(site_se[-ind])
       res<- rma(measure="GEN", yi=t(site_beta), sei=t(site_se),method="EE",control=list(maxiter=10000))
       effect_meta[i]<-exp(res$b)
       se_meta[i]<-exp(res$se)
       pvalue_meta[i]<-res$pval
       cilb[i]<-exp(res$ci.lb)
       ciub[i]<-exp(res$ci.ub)
       I2[i]<-res$I2
       het_p[i]<-res$QEp
       
       site_p<-as.numeric(site_p[-ind])
       p_count[i]<-length(which(site_p<=0.2)) 
    }
   else if (length(ind)==(nc-1)/4 | length(ind)==(nc-1)/4-1){ 
      effect_meta[i]<-NA
      pvalue_meta[i]<-NA
      se_meta[i]<-NA
      cilb[i]<-NA
      ciub[i]<-NA
      I2[i]<-NA
      het_p[i]<-NA
      p_count[i]<-NA
   }
 }
  re<-cbind(data,effect_meta,pvalue_meta,se_meta,cilb,ciub,I2,het_p,p_count)
  #write.csv(re, file = "07.Microbial_GWAS/test_mete_pvalue_cut15.csv",row.names=F)
  write.csv(re,dataout,row.names=F)
}


#######################################################################
##　多重检验校正


multi_adjust<-function(datain,dataout){

   #p_test<-read.csv("07.Microbial_GWAS/datasets/mel_vsv_os_meta_pvalue.csv",header=T)
   #p_test<-read.csv(datain,header=T)
   #data<-nsclc_resp
   data<-datain
   re_rnma<-subset(data,is.na(pvalue_meta)==F,drop=T)
   re_rnma<-re_rnma[order(re_rnma$X),]
   
   Spe<-unlist(strsplit(re_rnma$X[1],':'))[1]
   ad_p<-re_rnma$pvalue_meta[1]
   fdr_all<-0
   bonferroni_all<-0

   nr<-nrow(re_rnma)
   ###re_rnma<-re_rnma[order(re_rnma$X),]  

   #i<-1
  for(i in 2:nr){
     Species<-unlist(strsplit(re_rnma$X[i],':'))[1]
     if(Species==Spe){
        ad_p<-c(ad_p,re_rnma$pvalue_meta[i])
      }
     else if(Species!=Spe){
        fdr<-p.adjust(ad_p,method = 'fdr')
        bonferroni<-p.adjust(ad_p,method = 'bonferroni') 
        fdr_all<-c(fdr_all,fdr)
        bonferroni_all<-c(bonferroni_all,fdr)
        Spe<-unlist(strsplit(re_rnma$X[i],':'))[1]
        ad_p<-re_rnma$pvalue_meta[i]
       }
     if(i==nr){
       fdr<-p.adjust(ad_p,method = 'fdr')
       bonferroni<-p.adjust(ad_p,method = 'bonferroni') 
       fdr_all<-c(fdr_all,fdr)
       bonferroni_all<-c(bonferroni_all,fdr)
     }
  }

  re_out<-data.frame(cbind(re_rnma,fdr_all[-1],bonferroni_all[-1]))
  colnames(re_out)[c((ncol(re_out)-1),ncol(re_out))]<-c("fdr.p","bonf.p")
  #write.csv(re_out, file = "08.Microbial_GWAS/all_vsv_meta_pvalue_adjust.csv")
  write.csv(re_out, dataout,row.names=F)
}







average_matrix<-function(mat1,mat2){
   #mat1<-vsgv_msv_dist_i
   #mat2<-dsgv_msv_dist_i
   
   mat1_samples<-colnames(mat1)
   mat2_samples<-colnames(mat2)
  
   mat_samples<-union(mat1_samples,mat2_samples)
   nmat<-length(mat_samples)
   
   merge_mat<- matrix(NA,nrow=nmat,ncol=nmat,dimnames=list(mat_samples,mat_samples))

   for(i in 1:nmat){
       for(j in 1:nmat){
         if(mat_samples[i] %in% mat1_samples & mat_samples[j] %in% mat1_samples & mat_samples[i] %in% mat2_samples & mat_samples[j] %in% mat2_samples){
             merge_mat[mat_samples[i],mat_samples[j]]<-mean(c(mat1[mat_samples[i],mat_samples[j]],mat2[mat_samples[i],mat_samples[j]]),na.rm=T)
            } 
         else if (mat_samples[i] %in% mat1_samples & mat_samples[j] %in% mat1_samples & ((mat_samples[i] %in% mat2_samples)==F | (mat_samples[j] %in% mat2_samples)==F)){
             merge_mat[mat_samples[i],mat_samples[j]]<-mat1[mat_samples[i],mat_samples[j]]
            } 
        else if (mat_samples[i] %in% mat2_samples & mat_samples[j] %in% mat2_samples & ((mat_samples[i] %in% mat1_samples)==F | (mat_samples[j] %in% mat1_samples)==F)){
             merge_mat[mat_samples[i],mat_samples[j]]<-mat2[mat_samples[i],mat_samples[j]]
            }
         }
      }
   merge_mat
 }
 
 
 

dataset_msv_dist<-function(vsgv_in,dsgv_in,info_final,msv_dist_out,msv_dist_std_out){
  msv_dist<-NULL

  vsgv<-vsgv_in
  dsgv<-dsgv_in

  #vsgv<-vsgv_sub
  #dsgv<-dsgv_sub

  for (i in c(1:nrow(info_final))){
    #i<-16
    file_name<-str_replace_all(info_final$organism[i], "\\/", "\\.")
    spe_name<-str_replace_all(info_final$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    vsgv_i<-data.frame(vsgv[,grep(spe_name,colnames(vsgv))])
    dsgv_i<-data.frame(dsgv[,grep(spe_name,colnames(dsgv))])

    if(ncol(vsgv_i)>1 & ncol(vsgv_i)>1 ){
       vsgv_msv_i<-vsgv_i %>%
       na.omit(.)
       vsgv_msv_i = vsgv_msv_i[rowSums(vsgv_msv_i[])>0,]

       dsgv_msv_i<-dsgv_i %>%
       na.omit(.)
       dsgv_msv_i = dsgv_msv_i[rowSums(dsgv_msv_i[])>0,]

      if(nrow(vsgv_msv_i)>0 & nrow(dsgv_msv_i)>0){
        vsgv_msv_dist_i <- as.matrix(vegdist(as.data.frame(vsgv_msv_i),method = "canberra"))
        dsgv_msv_dist_i <- as.matrix(vegdist(as.data.frame(dsgv_msv_i),method = "jaccard"))
        msv_dist_i<-average_matrix(vsgv_msv_dist_i,dsgv_msv_dist_i)
        }
      else if(nrow(vsgv_msv_i)>0 & nrow(dsgv_msv_i)==0){
         vsgv_msv_dist_i <- as.matrix(vegdist(as.data.frame(vsgv_msv_i),method = "canberra"))
         msv_dist_i<-vsgv_msv_dist_i
         }
     else if(nrow(dsgv_msv_i)>0 & nrow(vsgv_msv_i)==0){
         dsgv_msv_dist_i <- as.matrix(vegdist(as.data.frame(dsgv_msv_i),method = "jaccard"))
         msv_dist_i<-dsgv_msv_dist_i
         }
      else {
         msv_dist_i<-NA
         }
       }
     else {
      msv_dist_i<-NA
    }
      msv_dist[[i]]<-msv_dist_i
   }

   names(msv_dist)<-paste('msv_',info_final$organism, sep = '')
   msv_dist_std <- lapply(msv_dist, myfun<-function(x){x/max(x,na.rm = T)})

   save(msv_dist, file = msv_dist_out)
   save(msv_dist_std, file = msv_dist_std_out)

   #save(msv_dist, file = "01.cleanData/SV_all/distMat/PRJEB22863_msv_dist.RData")
   #save(msv_dist_std, file = "01.cleanData/SV_all/distMat/PRJEB22863_msv_dist_std.RData")
 }





Q4_vsv_pfs<-function(datain,study){
  #pic_sub<-PetersBA_2020_pic
  pic_sub<-datain
  qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
    pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }
  ##return(pic_sub)
  
  ratio<-table(pic_sub$pfs_12_months,pic_sub$subtype)
  sum_vec<-apply(ratio,2,sum)
  ratio_new<-matrix(NA,2,4)

  for(i in 1:4){
   ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
  }
   colnames(ratio_new)<-c("Q1","Q2","Q3","Q4")
   rownames(ratio_new)<-c("PFS_less12","PFS_larger12")
   ratio_pic<-data.frame(ratio_new)

   ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
   ratio_pic<-melt(ratio_pic,id="chara")

   ratio_pic$study<-study
   return(ratio_pic)
}




Q4_vsv_resp<-function(datain,study){
  #pic_sub<-PetersBA_2020_pic
  pic_sub<-datain
  qt<-quantile(pic_sub$temp, probs = c(0.25,0.5,0.75),na.rm=T)
    pic_sub$subtype<-0
    for (j in 1:nrow(pic_sub)){
      if (is.na(pic_sub$temp[j])==T) {pic_sub$subtype[j]=NA}
      else if (pic_sub$temp[j]<=qt[1]){pic_sub$subtype[j]=1}
      else if (pic_sub$temp[j]>qt[1] & pic_sub$temp[j]<=qt[2]){pic_sub$subtype[j]=2}
      else if (pic_sub$temp[j]>qt[2] & pic_sub$temp[j]<=qt[3]){pic_sub$subtype[j]=3}
      else if (pic_sub$temp[j]>qt[3]){pic_sub$subtype[j]=4}
     }
  ##return(pic_sub)
  
  ratio<-table(pic_sub$response_code,pic_sub$subtype)
  sum_vec<-apply(ratio,2,sum)
  ratio_new<-matrix(NA,2,4)

  for(i in 1:4){
   ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
  }
   colnames(ratio_new)<-c("Q1","Q2","Q3","Q4")
   rownames(ratio_new)<-c("Non_RE","RE")
   ratio_pic<-data.frame(ratio_new)

   ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
   ratio_pic<-melt(ratio_pic,id="chara")

   ratio_pic$study<-study
   return(ratio_pic)
}


dsv_resp<-function(datain,study){

   ratio<-table(datain$response_code,datain$temp)
   sum_vec<-apply(ratio,2,sum)
   ratio_new<-matrix(NA,2,2)

  for(i in 1:2){
   ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
   }

  colnames(ratio_new)<-c("Delection","Non-delection")
  rownames(ratio_new)<-c("Non_RE","RE")

  #调整数据布局
  ratio_pic<-data.frame(t(ratio_new))
  ratio_pic$chara<-rownames(ratio_pic)

  ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
  ratio_pic<-melt(ratio_pic,id="chara")

  ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic
  #table(pic$temp)
  ratio_pic$study<-study
  
  return(ratio_pic)
}





dsv_pfs<-function(datain,study){
   
   #datain<-LeeKA_2022_pic
   ratio<-table(datain$pfs_12_months,datain$temp)
   sum_vec<-apply(ratio,2,sum)
   ratio_new<-matrix(NA,2,2)

  for(i in 1:2){
   ratio_new[,i]=round(ratio[,i]/sum_vec[i],2)
   }

  colnames(ratio_new)<-c("Delection","Non-delection")
  rownames(ratio_new)<-c("PFS_less12","PFS_larger12")

  #调整数据布局
  ratio_pic<-data.frame(t(ratio_new))
  ratio_pic$chara<-rownames(ratio_pic)

  ratio_pic$chara <- factor(rownames(ratio_pic), levels=rev(rownames(ratio_pic)))
  ratio_pic<-melt(ratio_pic,id="chara")

  ratio_pic %>% 
  group_by(chara) %>% 
  mutate(new_col=(1-cumsum(value))+0.05) -> ratio_pic
  #table(pic$temp)
  ratio_pic$study<-study
  
  return(ratio_pic)
}

