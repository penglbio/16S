####
library(grid) 
rm(list = ls())
wdir <- "/home/pengl/Documents/data_analysis/16s/2017_4_20/2017_4_3_rep/"
setwd(wdir)
source("R_script/tern_e.R")


inputFolder <- "./data_for_figure/"
outputFolder <- "./figure/"

#load design matrix and subset by experiment------------------------------------------------------
designfile <- "map.tsv"
map <- read.table(paste(inputFolder, designfile, sep = ""), sep = "\t",row.names = 1,header=T,stringsAsFactors = F) 


#ACM raw data  change to dat_norm_mat and dat_log_mat-------------------------------------------
# RA_value count function
RA_count<-function(m, threshold = 2, scale = 1000) {
  m <- m / rowSums(m, na=T) * scale
  return(m)
}
ACM<-read.table(paste(inputFolder,"ACM_norm_tax_form.txt",sep = ""),header=T,row.names = 1,sep="\t")
dim(ACM)
dat_norm_ACM<-t(RA_count(t(ACM[,1:84])))
dat_log_ACM <- log2(dat_norm_ACM+1)
###phyrum color
ACM_tax<-ACM[,85:91]
dim(ACM_tax)
ACM_tax[1:8,]
names(sort(table(droplevels(ACM_tax[,"Phylum"])),decr=T))
PHYLAcols <- as.character(ACM_tax[,"Phylum"])
names(PHYLAcols) <- rownames(ACM_tax)
PHYLAcols
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Proteobacteria"]] <- "seagreen4"     # "palegreen4" #"chartreuse4" "springgreen4" #
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Actinobacteria"]] <- "red"   #"tomato2"# "deeppink3"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Bacteroidetes"]] <- "steelblue2"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Acidobacteria"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Verrucomicrobia"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Firmicutes"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Gemmatimonadetes"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Chloroflexi"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Planctomycetes"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "TM7"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Nitrospirae"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "not_assigned"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Armatimonadetes"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Chlorobi"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "OD1"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "WS3"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Chlamydiae"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "FBP"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Fibrobacteres"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "Spirochaetes"]] <- "grey"
PHYLAcols[names(PHYLAcols)[PHYLAcols == "TM6"]] <- "grey"

##BS2 means and data
BS2<-rownames(map)[map$Species=="Fi_Bulk_soil"]
dat_RA_log_bs2<-dat_log_ACM[,BS2]
dat_RA_log_mean_bs2<-apply(dat_RA_log_bs2, 1, mean)
dat_RA_log_mean_bs2[dat_RA_log_mmean_bs2==0]<-0.0000001
##aov_analysis and p.adjust function
aov_adjust<-function(data){
  data1 <- data
  group<-as.factor(c(1,1,1,2,2,2))
  aov_file=NULL
  for(i in 1:length(data1[,1])){a=data1[i,];df<-data.frame(x=as.numeric(a),y=group);aov.out<-TukeyHSD(aov(x~y,data=df));pvalue<-aov.out$y[1,];aov_file<-rbind(aov_file,pvalue)}
  row.names(aov_file)<-row.names(data1)
  c_s_aov<-as.data.frame(aov_file)
  BH<-p.adjust(c_s_aov$`p adj`,"BH")
  c_s_aov<-cbind(c_s_aov,BH)
  return(c_s_aov)
}

##hen1-6 log2(BS2,rhizosphere,root) means and rootOTU
hen_id<-rownames(map)[map$Genotype=="hen1-6"]
dat_RA_log_hen<-dat_log_ACM[,hen_id]
dat_RA_log_mean_hen_rhizo<- apply(dat_RA_log_hen[,1:3], 1, mean)
dat_RA_log_mean_hen_rhizo[dat_RA_log_mean_hen_rhizo==0]<-0.0000001
dat_RA_log_mean_hen_root<- apply(dat_RA_log_hen[,4:6], 1, mean)
dat_RA_log_mean_hen_root[dat_RA_log_mean_hen_root==0]<-0.0000001
hen_RA_log_mean <- cbind(dat_RA_log_mean_hen_root,dat_RA_log_mean_hen_rhizo,dat_RA_log_mean_bs2)
colnames(hen_RA_log_mean)<-c("Root","Rhizosphere","BS2")
# hen_RA_log_mean[1:5,]
RA_log_hen_root<-dat_RA_log_hen[,4:6]
data <- cbind(dat_RA_log_bs2,RA_log_hen_root)
hen_aov<-aov_adjust(data)
hen_fold<-apply(dat_norm_ACM[,hen_id][,4:6],1,mean)/apply(dat_norm_ACM[,BS2],1,mean)
write.table(hen_aov,file = "hen_aov",sep="\t")
pid<-hen_fold[rownames(hen_aov)[which(hen_aov[,5]<=0.1)]]
hen_rootOTUs<-names(pid)[pid>=2|pid<=0.5]
names(sort(table(droplevels(ACM_tax[hen_rootOTUs,][,"Phylum"])),decr=T))
rownames(hen_RA_log_mean)
hen_rootOTUs_colors <- PHYLAcols[rownames(hen_RA_log_mean)]
hen_rootOTUs_colors[names(hen_rootOTUs_colors)[!names(hen_rootOTUs_colors) %in% hen_rootOTUs]] <- "grey"
pdf(paste(outputFolder,"rootOTUs_ternary_plots.pdf", sep=""), paper="special", width=5, height=5)
tern_e(hen_RA_log_mean, col=hen_rootOTUs_colors, prop=T, grid=T,  bg="transparent", pch=19, main="hen_rootOTUs",cex=0.5)
# dev.off()

##dcl2/3/4 log2(BS2,rhizosphere,root) means and rootOTU
dcl_234_id<-rownames(map)[map$Genotype=="dcl2/3/4"]
dcl_234_id
dat_RA_log_dcl_234<-dat_log_ACM[,dcl_234_id]
dat_RA_log_mean_dcl_234_rhizo<- apply(dat_RA_log_dcl_234[,1:3], 1, mean)
dat_RA_log_mean_dcl_234_rhizo[dat_RA_log_mean_dcl_234_rhizo==0]<-0.0000001
dat_RA_log_mean_dcl_234_root<- apply(dat_RA_log_dcl_234[,4:6], 1, mean)
dat_RA_log_mean_dcl_234_root[dat_RA_log_mean_dcl_234_root==0]<-0.0000001
dcl_234_RA_log_mean <- cbind(dat_RA_log_mean_dcl_234_root,dat_RA_log_mean_dcl_234_rhizo,dat_RA_log_mean_bs2)
colnames(dcl_234_RA_log_mean)<-c("Root","Rhizosphere","BS2")
# dcl_234_RA_log_mean[1:5,]
RA_log_dcl_234_root<-dat_RA_log_dcl_234[,4:6]
data <- cbind(dat_RA_log_bs2,RA_log_dcl_234_root)
dcl_234_aov<-aov_adjust(data)
dcl_234_fold<-apply(dat_norm_ACM[,dcl_234_id][,4:6],1,mean)/apply(dat_norm_ACM[,BS2],1,mean)
write.table(dcl_234_aov,file = "dcl_234_aov",sep="\t")
pid<-dcl_234_fold[rownames(dcl_234_aov)[which(dcl_234_aov[,5]<=0.1)]]
dcl_234_rootOTUs<-names(pid)[pid>=2|pid<=0.5]
names(sort(table(droplevels(ACM_tax[dcl_234_rootOTUs,][,"Phylum"])),decr=T))
rownames(dcl_234_RA_log_mean)
dcl_234_rootOTUs_colors <- PHYLAcols[rownames(dcl_234_RA_log_mean)]
dcl_234_rootOTUs_colors[names(dcl_234_rootOTUs_colors)[!names(dcl_234_rootOTUs_colors) %in% dcl_234_rootOTUs]] <- "grey"
# postscript(paste(outputFolder,"dcl_234_rootOTUs.pdf", sep=""), paper="special", width=5, height=5, horizontal = FALSE)
tern_e(dcl_234_RA_log_mean, col=dcl_234_rootOTUs_colors, prop=T, grid=T,  bg="transparent", pch=19, main="dcl_234_rootOTUs",cex=0.5)
# dev.off()


##dcl3 log2(BS2,rhizosphere,root) means and rootOTU
dcl3_id<-rownames(map)[map$Genotype=="dcl3"]
# dcl3_id
dat_RA_log_dcl3<-dat_log_ACM[,dcl3_id]
# dat_RA_log_dcl3[1:5,]
dat_RA_log_mean_dcl3_rhizo<- apply(dat_RA_log_dcl3[,1:3], 1, mean)
dat_RA_log_mean_dcl3_rhizo[dat_RA_log_mean_dcl3_rhizo==0]<-0.0000001
dat_RA_log_mean_dcl3_root<- apply(dat_RA_log_dcl3[,4:6], 1, mean)
dat_RA_log_mean_dcl3_root[dat_RA_log_mean_dcl3_root==0]<-0.0000001
dcl3_RA_log_mean <- cbind(dat_RA_log_mean_dcl3_root,dat_RA_log_mean_dcl3_rhizo,dat_RA_log_mean_bs2)
colnames(dcl3_RA_log_mean)<-c("Root","Rhizosphere","BS2")
# dcl3_RA_log_mean[1:5,]
RA_log_dcl3_root<-dat_RA_log_dcl3[,4:6]
data <- cbind(dat_RA_log_bs2,RA_log_dcl3_root)
dcl3_aov<-aov_adjust(data)
dcl3_fold<-apply(dat_norm_ACM[,dcl3_id][,4:6],1,mean)/apply(dat_norm_ACM[,BS2],1,mean)
write.table(dcl3_aov,file = "dcl3_aov",sep="\t")
pid<-dcl3_fold[rownames(dcl3_aov)[which(dcl3_aov[,5]<=0.1)]]
dcl3_rootOTUs<-names(pid)[pid>=2|pid<=0.5]
names(sort(table(droplevels(ACM_tax[dcl3_rootOTUs,][,"Phylum"])),decr=T))
rownames(dcl3_RA_log_mean)
dcl3_rootOTUs_colors <- PHYLAcols[rownames(dcl3_RA_log_mean)]
dcl3_rootOTUs_colors[names(dcl3_rootOTUs_colors)[!names(dcl3_rootOTUs_colors) %in% dcl3_rootOTUs]] <- "grey"
# postscript(paste(outputFolder,"dcl3_rootOTUs.pdf", sep=""), paper="special", width=5, height=5, horizontal = FALSE)
tern_e(dcl3_RA_log_mean, col=dcl3_rootOTUs_colors, prop=T, grid=T,  bg="transparent", pch=19, main="dcl3_rootOTUs",cex=0.5)
dev.off()
