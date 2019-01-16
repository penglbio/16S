####this program is for analysis the Fold change between the samples


rm(list = ls())
wdir <- "/home/pengl/Documents/data_analysis/16s/2017_4_20/2017_4_3_rep/"
setwd(wdir)
inputFolder <- "./data_for_stas_analysis/"
outputFolder <- "./stats_result/one_by_one_aov/"
library("dplyr")
library("tidyr")
##expermant design
map<-read.table(paste(inputFolder,"map.tsv",sep=""),header=T,row.names=1,sep = "\t",stringsAsFactors = F)
BS1<-"In_Bulk_soil"
BS2<-"Fi_Bulk_soil"
rhizo_wt<-"rhizosphere_wt"
root_wt<-"root_wt"
rhizo_mutant<-unique(map$Species[grep("rhizosphere_mutant*",map$Species)])
root_mutant<-unique(map$Species[grep("root_mutant*",map$Species)])
rhizo_root<-c(rhizo_wt,rhizo_mutant,root_wt,root_mutant)

##FC
FC_count<-function(control_name,sample_name,data){
  control_id<-rownames(map)[map$Species==control_name]
  sample_id<-rownames(map)[map$Species==sample_name]
  data<-data
  data1<-data[,c(control_id,sample_id)]
  group<-c(1,1,1,2,2,2)
  apply(data1,1,function(x){otu<-data.frame(x=as.numeric(x),y=group);mean_c_s<-aggregate(otu[,1],otu[,2],mean);FC <- mean_c_s[2]/mean_c_s[1]})
}

##logFC---------------------------------------------------------------------------------------------
#wt_vs_mutant_rhizo
pwm<-read.delim(paste(inputFolder,"wt_vs_mutant_rhizo.txt", sep=""),header=T,row.names=1,sep = "\t")
dim(pwm)
data<-pwm[,3:41]
comp_sam_FC<-data_frame()
  for (i  in rhizo_mutant ) {
  FC_line<-FC_count(rhizo_wt,i,data )
  FC_line<-as.data.frame(FC_line)
  if(nrow(comp_sam_FC)==0)
  {comp_sam_FC<-FC_line}
  else
  {comp_sam_FC<-cbind(comp_sam_FC,FC_line)}
  }
map$Genotype[map$Species %in% rhizo_mutant]
colnames(comp_sam_FC)<-c(unique(map$Genotype[map$Species %in% rhizo_mutant]))
otu_rhizo_wt_mutant_FC<- comp_sam_FC

##wt_vs_mutant_root
pwm<-read.delim(paste(inputFolder,"wt_vs_mutant_root.txt", sep=""),header=T,row.names=1,sep = "\t")
dim(pwm)
data<-pwm[,2:40]
comp_sam_FC<-data_frame()
for (i  in root_mutant ) {
  FC_line<-FC_count(root_wt,i,data )
  FC_line<-as.data.frame(FC_line)
  if(nrow(comp_sam_FC)==0)
  {comp_sam_FC<-FC_line}
  else
  {comp_sam_FC<-cbind(comp_sam_FC,FC_line)}
}
colnames(comp_sam_FC)<-c(unique(map$Genotype[map$Species %in% root_mutant]))
otu_root_wt_mutant_FC<- comp_sam_FC
##bs1_bs2
pwm<-read.delim(paste(inputFolder,"BS1_vs_BS2.txt", sep=""),header=T,row.names=1,sep = "\t")
dim(pwm)
data<-pwm[,2:7]
comp_sam_FC<-data_frame()
for (i  in BS2 ) {
  FC_line<-FC_count(BS1,i,data )
  FC_line<-as.data.frame(FC_line)
  if(nrow(comp_sam_FC)==0)
  {comp_sam_FC<-FC_line}
  else
  {comp_sam_FC<-cbind(comp_sam_FC,FC_line)}
}
colnames(comp_sam_FC)<-c(unique(map$Genotype[map$Species %in% BS2]))
BS1_vs_BS2_FC<- comp_sam_FC



BS1_id<-rownames(map)[grep("In_Bulk_soil", map$Species)]
BS2_id<-rownames(map)[grep("Fi_Bulk_soil", map$Species)]
rhizo_wt_id<-rownames(map)[grep("rhizosphere_wt", map$Species)]
root_wt_id<-rownames(map)[grep("root_wt", map$Species)]
rhizo_mutant_id<-rownames(map)[grep("rhizosphere_mutant*", map$Species)]
root_mutant_id<-rownames(map)[grep("root_mutant*", map$Species)]

##phyrum rhizo_wt_mutant
pwm<-read.delim(paste(inputFolder,"wt_vs_rhizo_phyrum.txt", sep=""),header=F,row.names=1,sep = "\t")
colnames(pwm)<-c(rhizo_wt_id,rhizo_mutant_id)
data<-pwm
comp_sam_FC<-data_frame()
for (i  in rhizo_mutant ) {
  FC_line<-FC_count(rhizo_wt,i,data )
  FC_line<-as.data.frame(FC_line)
  if(nrow(comp_sam_FC)==0)
  {comp_sam_FC<-FC_line}
  else
  {comp_sam_FC<-cbind(comp_sam_FC,FC_line)}
}
colnames(comp_sam_FC)<-c(unique(map$Genotype[map$Species %in% rhizo_mutant]))
phyrum_rhizo_wt_mutant_FC<- comp_sam_FC

##phyrum root_wt_mutant
pwm<-read.delim(paste(inputFolder,"wt_vs_root_phyrum.txt", sep=""),header=F,row.names=1,sep = "\t")
colnames(pwm)<-c(root_wt_id,root_mutant_id)
data<-pwm
comp_sam_FC<-data_frame()
for (i  in root_mutant ) {
  FC_line<-FC_count(root_wt,i,data )
  FC_line<-as.data.frame(FC_line)
  if(nrow(comp_sam_FC)==0)
  {comp_sam_FC<-FC_line}
  else
  {comp_sam_FC<-cbind(comp_sam_FC,FC_line)}
}
colnames(comp_sam_FC)<-c(unique(map$Genotype[map$Species %in% root_mutant]))
phyrum_root_wt_mutant_FC<- comp_sam_FC


##phyrum BS2_rhizo_root
pwm<-read.delim(paste(inputFolder,"BS_vs_rhizo_and_root_phyrum.txt", sep=""),header=F,row.names=1,sep = "\t")
colnames(pwm)<-c(BS2_id,rhizo_wt_id,rhizo_mutant_id,root_wt_id,root_mutant_id)
data<-pwm
comp_sam_FC<-data_frame()
for (i  in rhizo_root ) {
  FC_line<-FC_count(BS2,i,data )
  FC_line<-as.data.frame(FC_line)
  if(nrow(comp_sam_FC)==0)
  {comp_sam_FC<-FC_line}
  else
  {comp_sam_FC<-cbind(comp_sam_FC,FC_line)}
}
colnames(comp_sam_FC)<-c(unique(map$Genotype[map$Species %in% rhizo_root]))
phyrum_rhizo_root_FC<- comp_sam_FC





