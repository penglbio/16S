argv <- commandArgs(T)

sample_stat<-argv[1]
barcode_seq<-argv[2]

data<-read.table(sample_stat)

data$percentage=100*data$V2/sum(data$V2)

colnames(data)<-c("seq","count","precentage")

sample<-read.table(barcode_seq)

colnames(sample) <- c("sample","seq")

sample_data <- merge(sample,data,by.x=2,by.y=1,all.x=T)

sample_data<-sample_data[order(sample_data$sample),]

write.table(sample_data,file=sample_stat,quote=F,sep="\t",row.names=F)
