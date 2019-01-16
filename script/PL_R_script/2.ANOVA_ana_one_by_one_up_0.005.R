###ANOVA ana one by one RA up 0.005
rm(list = ls())
wdir <- "/var/data/19.Richa_new/"
setwd(wdir)
inputFolder <- "./3.data_for_figure/"
outputFolder <- "./2.stastics_analysis_result/"
source("script/Cbind.R")
library("dplyr")
library("tidyr")
library("plyr")
library("qpcR")
##aov_analysis and p.adjust function
aov_adjust<-function(control_name,sample_name,data){
  control_id<-rownames(map)[map$Genotype==control_name]
  sample_id<-rownames(map)[map$Genotype==sample_name]
  data[,c(control_id,sample_id)]
  data1<-data[,c(control_id,sample_id)]
  data1<-log2(data1+1)
  dim(data1)
  data[,2]
  group<-as.factor( c(rep(1,length(control_id)),rep(2,length(sample_id))))
  aov_file=NULL
  for(i in 1:length(data1[,1])){a=data1[i,];df<-data.frame(x=as.numeric(a),y=group);aov.out<-TukeyHSD(aov(x~y,data=df));pvalue<-aov.out$y[1,];aov_file<-rbind(aov_file,pvalue)}
  row.names(aov_file)<-row.names(data1)
  c_s_aov<-as.data.frame(aov_file)
  BH<-p.adjust(c_s_aov$`p adj`,"fdr")
  c_s_aov<-cbind(c_s_aov,BH=BH)
  fold<-apply(data[,sample_id],1,mean)/apply(data[,control_id],1,mean)
  c_s_aov<-cbind(c_s_aov,FC=fold)
  return(c_s_aov)
}


cmp_wt_mutant <- function(aov_file,s_control,s_mutant,data,type){
  aov_file<-data.frame()
  for (i in s_mutant)
  {
    aov_df<-aov_adjust(s_control,i,data)
    if( nrow(aov_file)==0)
    {aov_file<-aov_df}
    else
    {aov_file<-cbind(aov_file,aov_df)}
  }
  colnames(aov_file)<-paste(rep(s_mutant,each=6),colnames(aov_file),sep="_")
  if(type=="otu"){
    aov_file<-data.frame(aov_file,tax)
  }
  else{
    rownames(aov_file)<-rownames(data)
    aov_file
  }
}

##read_tables(acm phylum family)
##read_table
argv <- commandArgs(T)
otu_table<-"otu_change_syringae.txt"
ACM_table<-"ACM_norm_tax_form_change_syringae.txt"
ACM_RA<-"ACM_RA_norm_tax_form.txt"
fam_RA<-"Family_RA_table.txt"
phy_RA<-"Phy_RA_table.txt"

otu_table<-read.delim(paste(inputFolder,otu_table, sep=""),header=T,row.names= 1,sep = "\t",stringsAsFactors = F)
acm_table<-read.table(paste(inputFolder,ACM_table, sep=""),header=T,row.names = 1,sep = "\t",stringsAsFactors = F)
#group_type <- argv[6]
##expermant design

map<-read.table(paste(inputFolder,"map.tsv",sep=""),header=T,row.names=1,sep = "\t",stringsAsFactors = F)
rownames(map) <- paste0("X",rownames(map))
map <- map[rownames(map[rownames(map)%in%colnames(otu_table),]),]
Rhizo <- rownames(map[grepl("rizho",map$Genotype),])
root <- rownames(map[grepl("root",map$Genotype),])
soil <- rownames(map[grepl("soil",map$Genotype),])

rhizo_wt <- unique(map[Rhizo,]$Genotype)[1]
rhizo_mutant<- unique(map[Rhizo,]$Genotype)[-1]

root_wt <- unique(map[root,]$Genotype)[1]
root_mutant<- unique(map[root,]$Genotype)[-1]

bs1 <- unique(map[soil,]$Genotype)[1]
bs2<- unique(map[soil,]$Genotype)[-1]



#--------ACM------------
pwm<-read.delim(paste(inputFolder,ACM_RA, sep=""),header=T,row.names=1,sep = "\t",stringsAsFactors = F)
data <- pwm[,-((ncol(pwm)-6):ncol(pwm))]
data_5<-apply(data,2,function(x){ifelse(x>=5,1,0)})
ACM_RA_5<-data[apply(data_5,1,sum)>=1,]
tax<-apply(pwm[rownames(ACM_RA_5),((ncol(pwm)-6):ncol(pwm))],1,paste0,collapse=";")

#--------Family-----------
fam_pwm<-read.delim(paste(inputFolder,fam_RA, sep=""),header=T,row.names=1,sep = "\t")
fam_data<-fam_pwm[,2:ncol(fam_pwm)]
fam_data_5<-apply(fam_data,2,function(x){ifelse(x>=5,1,0)})
fam_RA_5<-fam_data[apply(fam_data_5,1,sum)>=1,]

#-------phylum-----------
phy_pwm<-read.delim(paste(inputFolder,phy_RA, sep=""),header=T,row.names=1,sep = "\t")
phy_data<-phy_pwm
phy_data_5<-apply(phy_data,2,function(x){ifelse(x>=5,1,0)})
phy_RA_5<-phy_data[apply(phy_data_5,1,sum)>=1,]

#otu_wt_mutant
otu_rhizo_wt_mutant <- cmp_wt_mutant(rhizo_wt_mutant,rhizo_wt,rhizo_mutant,ACM_RA_5,"otu")
otu_root_wt_mutant <- cmp_wt_mutant(root_wt_mutant,root_wt,root_mutant,ACM_RA_5,"otu")
otu_bs1_bs2 <-cmp_wt_mutant(bs_wt_mutant,bs1,bs2,ACM_RA_5,"otu")

#family_wt_mutant
fam_rhizo_wt_mutant <- cmp_wt_mutant(rhizo_wt_mutant,rhizo_wt,rhizo_mutant,fam_RA_5,"fam")
fam_root_wt_mutant <- cmp_wt_mutant(root_wt_mutant,root_wt,root_mutant,fam_RA_5,"fam")
fam_bs1_bs2 <-cmp_wt_mutant(bs_wt_mutant,bs1,bs2,fam_RA_5,"fam")

##phylum_wt_mutant
phy_rhizo_wt_mutant <- cmp_wt_mutant(rhizo_wt_mutant,rhizo_wt,rhizo_mutant,phy_RA_5,"phy")
phy_root_wt_mutant <- cmp_wt_mutant(root_wt_mutant,root_wt,root_mutant,phy_RA_5,"phy")
phy_bs1_bs2 <-cmp_wt_mutant(bs_wt_mutant,bs1,bs2,phy_RA_5,"phy")

###ACM Fam phy RA>5
ACM_RA_5_with_Desc<-rbind(OTU_id=as.vector(map[colnames(ACM_RA_5),]$Ge),ACM_RA_5)
fam_RA_5_Desc<-rbind(OTU_id=as.vector(map[colnames(fam_RA_5),]$Genotype),fam_RA_5)
phy_RA_5_Desc<-rbind(OTU_id=as.vector(map[colnames(phy_RA_5),]$Genotype),phy_RA_5)

##write to excel
jgc <- function()
{
  gc()
  .jcall("java/lang/System", method = "gc")
}  

library("rJava")
library("xlsx")
set.seed(19790801)
n_sheets <- 40
options(java.parameters = "-Xmx8000m")
#df<-list(otu_table,acm_table,pwm,fam_pwm,phy_pwm,ACM_RA_5_with_geno,Fam_ACM_RA_5_geno,phy_ACM_RA_5_geno,rhizo_wt_mutant_with_tax,silique_wt_mutant_with_tax,root_wt_mutant_with_tax,stem_wt_mutant_with_tax,bs1_bs2_tax,
#family_rhizo_wt_mutant,family_silique_wt_mutant,family_root_wt_mutant,family_stem_wt_mutant,family_BS2_BS1,phy_rhizo_wt_mutant,phy_silique_wt_mutant,phy_root_wt_mutant,phy_stem_wt_mutant,phy_BS2_BS1)
wb<-createWorkbook()
sname <- c("ACM_table","ACM_RA","fam_RA","phy_RA","ACM_RA_5","fam_RA_5","phy_RA_5","otu_rhizo_wt_mutant_with_tax","otu_root_wt_mutant_with_tax",
          "otu_bs1_bs2_tax","family_rhizo_wt_mutant","family_root_wt_mutant","family_bs2_bs1",
           "phy_rhizo_wt_mutant","phy_root_wt_mutant","phy_bs2_bs1")
df<-list(acm_table,pwm,fam_pwm,phy_pwm,ACM_RA_5_with_Desc,fam_RA_5_Desc,phy_RA_5_Desc,otu_rhizo_wt_mutant,otu_root_wt_mutant,otu_bs1_bs2,
         fam_rhizo_wt_mutant,fam_root_wt_mutant,fam_bs1_bs2,
         phy_rhizo_wt_mutant,phy_root_wt_mutant,phy_bs1_bs2)
for (i in 1:16){
  jgc()
  sheet<-createSheet(wb,sheetName = as.character(sname[i]) )
  addDataFrame(df[[i]],sheet,row.names = TRUE)
}
saveWorkbook(wb,paste(outputFolder,"stat_otus_family_phyrum_result_with_FC_and_RA_table_group1.xlsx",sep=""))



