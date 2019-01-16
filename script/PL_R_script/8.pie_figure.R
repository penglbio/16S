rm(list=ls())
wdir <- "/home/pengl/Documents/data_analysis/16s/2017_4_20/2017_4_3_rep/"
setwd(wdir)
inputFolder <- "./data_for_figure/"
outputFolder <- "./figure/"
pie_data<-read.table(paste(inputFolder,"Pie_figure.txt",sep=""),header=T,sep = "\t",stringsAsFactors = F)
order<-table(pie_data$Order)
a<-unique(pie_data$Order)
order<-order[a]
fam<-table(pie_data$Family)
b<-unique(pie_data$Family)
fam<-fam[b]
gen<-table(pie_data$Genus)
gen["not assigned"]<-1
a<-2
names(a)<-"not assigned"
gen1<-c(gen[1:7],a,gen[8:15])
gen1
par(fig=c(0,1,0,1),cex=0.6)
pie(gen1,border = 0,radius = 0.9,init.angle = 200,col=c("indianred3","seagreen4","seagreen4",rep("steelblue2",13)))
par(fig=c(0,1,0,1),new=TRUE,cex=0.7)
fam
pie(fam,border = 0,radius = 0.4,init.angle = 200,col =c("red3","green4","green4",rep("blue2",10)))
legend("topright",legend=c("Actinobacteria","Bacteroidetes","Proteobacteria"),col=c("indianred3","seagreen4","steelblue2"),pch=16,bty = "n")
