argv <- commandArgs(T)
otu_table<-argv[1]
ACM_table<-argv[2]


data<-read.table(otu_table,header=T,row.names="otu_id",sep="\t")
data_wo_tax<-data[,-ncol(data)]
ACM<-apply(data_wo_tax,2,function(x){ifelse(x>=20,1,0)})
ACM_20<-data[apply(ACM,1,sum)>=1,]
write.table(ACM_20,ACM_table,sep="\t",quote=FALSE)
