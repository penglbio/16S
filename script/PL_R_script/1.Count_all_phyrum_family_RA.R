#ANOVA ANALYSIS for family
rm(list = ls())
getwd()
argv <- commandArgs(T)
argv[1]<-"/var/data/21.Richa_16s/"
argv[2]<-"otus_table_flt_sample.txt"
argv[3]<-"ACM_group1_norm_tax_form.txt"
argv[4]<-"Beta_diversity_TIC_w_unifrac.pdf"



wdir <- argv[1]
setwd(wdir)
input_file1 <- argv[2]
input_file2 <- argv[3]
outfig <- argv[4]



inputFolder <- "./1.data_for_statstics/"
outputFolder <- "./3.data_for_figure/"

source("./script/Cbind.R")
library("dplyr")
library("tidyr")
library("plyr")
library("qpcR")
##expermant design
map<-read.table(paste(inputFolder,"map.tsv",sep=""),header=T,row.names=1,sep = "\t")
rownames(map) <- paste0("X",rownames(map))
###change OTU talbe syringae
otu<-read.delim(paste(inputFolder,input_file1, sep=""),header=T,row.names=1,sep ="\t")
map <- map[rownames(map[rownames(map)%in%colnames(otu),]),]
otu1<-otu[grepl("syringae",otu$taxonomy),]
otu2<-otu[!grepl("syringae",otu$taxonomy),]
otu1$taxonomy<-sub("Enterobacteriales","Pseudomonadales",otu1$taxonomy)
otu1$taxonomy<-sub("Enterobacteriaceae","Pseudomonadaceae",otu1$taxonomy)
otu<-rbind(otu1,otu2)
write.table(otu,paste(outputFolder,"otu_change_syringae.txt",sep=""),row.names = T,col.names = T,sep="\t")

###count RA talbe
pwm<-read.delim(paste(inputFolder,input_file2, sep=""),header=T,row.names=1,sep ="\t")
dim(pwm)
pwm_c1<-pwm[pwm$Species=="syringae",]
pwm_c2<-pwm[pwm$Species!="syringae",]
pwm_c1$Order<-sub("Enterobacteriales","Pseudomonadales",pwm_c1$Order)
pwm_c1$Family<-sub("Enterobacteriaceae","Pseudomonadaceae",pwm_c1$Family)
pwm<-rbind(pwm_c1,pwm_c2)
ACM_RA<-apply(pwm[,-((ncol(pwm)-6):ncol(pwm))],2,function(x){as.numeric(x)/sum(as.numeric(x))*1000})
ACM_RA<-as.data.frame(ACM_RA)
##ACM_RA_with_tax_form.txt
tax<-pwm[,((ncol(pwm)-6):ncol(pwm))]
ACM_RA_with_tax<-cbind(ACM_RA,tax)
write.table(ACM_RA_with_tax,paste(outputFolder,"ACM_RA_norm_tax_form.txt",sep=""),row.names = T,col.names = T,sep="\t",quote = FALSE)
write.table(pwm,paste(outputFolder,"ACM_norm_tax_form_change_syringae.txt",sep=""),row.names = T,col.names = T,sep="\t",quote = FALSE)

#####count_family_RA_table
ACM_RA_with_tax<-cbind(ACM_RA,phyrum=pwm[,(ncol(pwm)-5)],Family=pwm[,(ncol(pwm)-2)])
pwm_assi_Fam<-ACM_RA_with_tax[which(!grepl("not_assigned", ACM_RA_with_tax$Family)),]
pwm_noassi_Fam<-ACM_RA_with_tax[which(grepl("not_assigned", ACM_RA_with_tax$Family)),]
Fam_RA<-rowsum(pwm_assi_Fam[,-((ncol(pwm_assi_Fam)-1):ncol(pwm_assi_Fam))],group =as.factor(paste(pwm_assi_Fam$phyrum,pwm_assi_Fam$Family,sep="__")))
Fam_RA_noassi<-rowsum(pwm_noassi_Fam[,-((ncol(pwm_assi_Fam)-1):ncol(pwm_assi_Fam))],group = as.factor(paste(pwm_noassi_Fam$Family,pwm_noassi_Fam$Family,sep="__")))
Fam_RA1<-rbind(Fam_RA,Fam_RA_noassi)
class(Fam_RA1)
Phy=NULL
Fam=NULL
for(i in 1:nrow(Fam_RA1)){Phy<-c(Phy,strsplit(rownames(Fam_RA1),"__")[[i]][1])}
for(i in 1:nrow(Fam_RA1)){Fam<-c(Fam,strsplit(rownames(Fam_RA1),"__")[[i]][2])}
Fam_RA2<-cbind(Fam_RA1,Phyrum=Phy,Family=Fam)
Fam_RA2<-data.frame(Family=Fam_RA2[,ncol(Fam_RA2)],Fam_RA2[,-((ncol(Fam_RA2)-1):ncol(Fam_RA2))])
write.table(Fam_RA2,paste(outputFolder,"Family_RA_table.txt",sep=""),row.names = T,col.names = T,sep="\t",quote = FALSE)

###count_phyrum_RA_table
ACM_RA_with_tax<-cbind(ACM_RA,phyrum=pwm[,(ncol(pwm)-5)])
Phy_RA<-rowsum(ACM_RA_with_tax[,-((ncol(pwm)-6):ncol(pwm))],group =as.factor(ACM_RA_with_tax$phyrum))
write.table(Phy_RA,paste(outputFolder,"Phy_RA_table.txt",sep=""),row.names = T,col.names = T,sep="\t",quote = FALSE)

