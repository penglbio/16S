rm(list = ls())
library(readxl)
library(RColorBrewer)

argv <- commandArgs(T)
wdir<-argv[1]
rm_sample<-as.numeric(strsplit(argv[2],"-")[[1]])
stat_otu_table<-argv[3]

sheet_name <- argv[4]
fig_name<-argv[5]
sam1 <- argv[6]
sam2 <- argv[7]
threads<-as.numeric(strsplit(argv[8][[1]],",")[[1]])

setwd(wdir)
source("script/PL_R_script/5.cpcoa_fun.R")
inputFolder <- "./3.data_for_figure/"
outputFolder <- "./4.figures/"
map<-read.table(paste(inputFolder,"map.tsv",sep=""),header=T,row.names=1,sep = "\t",stringsAsFactors = F)
map <- map[-c(rm_sample),]
rownames(map) <- paste0("X",rownames(map))
Sample <- map[grepl(sam1 ,map$Description),]
Sample <- rownames(Sample[grepl(sam2,Sample$Description),])

ACM_RA<-read.table(paste(inputFolder,"ACM_RA_norm_tax_form.txt",sep = ""),header=T,row.names = 1,sep="\t",stringsAsFactors = F)
RA_dat <- ACM_RA[,1:(ncol(ACM_RA)-7)]
data_5<-apply(RA_dat,2,function(x){ifelse(x>=5,1,0)})
RA_dat<-RA_dat[apply(data_5,1,sum)>=1,]
dat<-log2(RA_dat+1)
tax_20 <- ACM_RA[,((ncol(ACM_RA)-6):ncol(ACM_RA))]
##select sig_OTU rhizoOTU


sig_otu_select<- function(sig_sheet){
  # sig_sheet <- sheet_name
  data <- as.data.frame(read_excel(stat_otu_table,sheet = sig_sheet))
  rownames(data) <- data[,1]
  BH <- data[,grepl("BH",colnames(data))]
  FC <- data[,grepl("FC",colnames(data))]
  if(class(BH)=="data.frame")
  {
    rza <- apply(BH,2,function(x){ifelse(x<=threads[1],1,0)})
    rza[is.na(rza)]=0
    rzb <- apply(FC,2,function(x){ifelse(x<=threads[2]|x>=threads[3],1,0)})
    rzb[is.na(rzb)]=0
    rza1 <- rownames(data[apply(rza,1,sum)>=1,])
    rzb1 <- rownames(data[apply(rzb,1,sum)>=1,])
    OTU <- intersect(rza1,rzb1)
    fig_phy <- names(sort(table(tax_20[OTU,]$Phylum),decreasing = T))
    fig_phy <- fig_phy[grepl("[^not_assigned]",fig_phy)]
    if(length(fig_phy)>=5){
      fig_phy <- fig_phy[1:5]
    }
    
    OTU<-OTU[tax_20[OTU,]$Phylum %in% fig_phy]
    OTU_fig_phy <- list()
    OTU_fig_phy <- list(OTU,fig_phy) 
  }
  else{
    rza <- ifelse(BH<=threads[1],1,0)
    rza[is.na(rza)]=0
    rzb <- ifelse(FC<=threads[2]|FC>=threads[3],1,0)
    rzb[is.na(rzb)]=0
    names(rza) <- rownames(data)
    names(rzb) <- rownames(data)
    rza1 <- rownames(data[rza>=1,])
    rzb1 <- rownames(data[rzb>=1,])
    OTU <- intersect(rza1,rzb1)
    fig_phy <- names(sort(table(tax_20[OTU,]$Phylum),decreasing = T))
    fig_phy <- fig_phy[grepl("[^not_assigned]",fig_phy)]
    if(length(fig_phy)>=5){
      fig_phy <- fig_phy[1:5]
    }
    OTU<-OTU[tax_20[OTU,]$Phylum %in% fig_phy]
    OTU_fig_phy <- list()
    OTU_fig_phy <- list(OTU,fig_phy) 
  }
}
# extract abundance matrices by experiment
pcoa_bioplot<- function(Sample,sig_OTU,fig_phy,fig_name){
  # Sample <- Sample
  # sig_OTU <- sig_OTU_phy[[1]]
  # fig_phy <- sig_OTU_phy[[2]]
  # fig_name <- "Rhizo_sig.pdf"
  d_sample<-as.data.frame(map[Sample ,]$Description)
  rownames(d_sample)<-Sample
  colnames(d_sample)<-"Species"
  dat_OTU<-t(dat[sig_OTU,][,Sample])
  # perform the CAP analysis using the capscale function within the vegan package
  sqrt_transform = F
  cap <- capscale(dat_OTU ~ Species, data = d_sample, add = F,sqrt.dis = sqrt_transform)
  # perform the ANOVA-based permutation tests for each transformation
  perm_anova <- anova.cca(cap)
  # generate variance tables
  var_tbl <- variability_table(cap)
  cap_var<- cap_var_props(cap)
  # calculate confidence intervals for the variance of each constrained transf.
  ci <-  cca_ci(cap)
  # obtain the coordinates of the samples for each transformation (sample scores)
  wa <- cap$CCA$wa[, 1:2]
  # obtain the coordinates of the OTUs for each transformation (species scores)
  v <- cap$CCA$v[, 1:2]
  # get the centroids of the constrained factors (spp.)
  centroids <- cap$CCA$centroids[, 1:2]
  centroids <- centroids[paste0("Species",unique(d_sample$Species)),]
  # plot the different CPCoA transformations
  par(mar=c(2,4,3,10),xpd=TRUE)
  library(RColorBrewer)
  ful_pch <- 15:25
  ful_col <-brewer.pal(9,"Set1")
  plot_cap(p = wa, d = d_sample, col_var = unique(d_sample$Species), 
           pch_var = unique(d_sample$Species), 
           col_comp = rep("Species",length(d_sample$Species)), 
           pch_comp = rep("Species",length(d_sample$Species)), 
           shapes = ful_pch[1:length(unique(d_sample$Species))], 
           colors = ful_col[1:length(unique(d_sample$Species))], 
           file_name = paste(outputFolder, "PCoA_",fig_name, sep = ""),
           svg = F, constraint = "Species", ci = ci,
           var_tbl = var_tbl, cap_var = cap_var,
           perm_anova = perm_anova, cex.main = 0.8, cex = 1,mar_var=c(6,6,5,6),leg_x=0.32,leg_y=0.05,leg_point_type=ful_pch[1:length(unique(d_sample$Species))],
           leg_var=sub(".+?_","",rownames(centroids)))
  # plot biplots and species scores
  #labels_var should be set by the centroids.root rownames turn
  plot_biplots(p = v, centroids = centroids, pch=16, col="black",
               file_name = paste(outputFolder, "bioplot_",fig_name, sep=""), 
               svg = F, constraint="Species", ci = ci, 
               var_tbl = var_tbl, cap_var = cap_var, 
               perm_anova = perm_anova, abundances = dat_OTU, phy=fig_phy,
               otu_subset = NULL, taxonomy = tax_20, cex.main = 0.8, cex = 1,
               labels_var=sub(".+?_","",rownames(centroids)),
               mar_var = c(6,6,5,8),leg_x = 0.25,leg_y = 0.3,leg_var = fig_phy)
}

sig_OTU_phy <- sig_otu_select(sheet_name)
pcoa_bioplot(Sample,sig_OTU_phy[[1]],sig_OTU_phy[[2]],fig_name)

