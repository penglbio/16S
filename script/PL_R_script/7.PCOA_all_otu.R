rm(list = ls())
library(readxl)
library(RColorBrewer)
source("script/PL_R_script/7.cap_fun.R")
argv <- commandArgs(T)
   # argv[1] <- "/var/data/03_ZHM/16s/shenglan/" 
   # argv[2] <-  "20-25-30" 
   # argv[3] <- "root.pdf"
   # argv[4] <- "Root"
   # argv[5] <- "gp1|gp1_3FLAG|35S|Col"
wdir<-argv[1]
rm_sample <- as.numeric(strsplit(argv[2],"-")[[1]])
fig_name<-argv[3]
sam1 <- argv[4]
sam2 <- argv[5]

setwd(wdir)
inputFolder <- "3.data_for_figure/"
outputFolder <- "4.figures/"
ACM_RA<-read.table(paste(inputFolder,"ACM_RA_norm_tax_form.txt",sep = ""),header=T,row.names = 1,sep="\t")
#load design matrix and subset by experiment
map <- read.table(paste(inputFolder, "map.tsv", sep = ""),sep = "\t",row.names = 1,header=T,stringsAsFactors = F) 
map <- map[-c(rm_sample),]
rownames(map) <- paste0("X",rownames(map))
Sample <- map[grepl(sam1,map$Description),]
Sample <- rownames(Sample[grepl(sam2,Sample$Description),])

# load OTU abundance matrix 
RA_dat<-ACM_RA[,grep("X",colnames(ACM_RA))]
data_5<-apply(RA_dat,2,function(x){ifelse(x>=5,1,0)})
RA_dat<-RA_dat[apply(data_5,1,sum)>=1,]
dat<-log2(RA_dat+1)
##matrix select 



# perform the CAP analysis using the capscale function within the vegan package
pcoa<- function(Sample,fig_name){
  Sample <- Sample
  dat_OTU<-t(dat[,Sample])
  d_sample<-data.frame(Species=map[Sample,]$Description,row.names = Sample,stringsAsFactors = F)
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
  ful_pch <- 15:25
  ful_col <-brewer.pal(9,"Set1")
  plot_cap(p = wa, d = d_sample,  centroids = centroids,
           col_var = unique(d_sample$Species), 
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
}
pcoa(Sample,fig_name)
