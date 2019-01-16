data<-read.table("8.rarefaction_analysis/OTU_table_TIC/rarefaction_29000.txt",header=T,row.names="OTU_ID",sep="\t")
data_wo_tax<-data[,-79]
ACM<-apply(data_wo_tax,2,function(x){ifelse(x>=20,1,0)})
ACM_20<-data[apply(ACM,1,sum)>=1,]
write.table(ACM_20,"8.rarefaction_analysis/OTU_table_ACM/ACM_table.txt",,sep="\t")
