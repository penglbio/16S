# Libraries and functions
library("ShotgunFunctionalizeR")    #  from http://shotgun.math.chalmers.se
library("vegan")
source("script/PL_R_script/4.emiel_diversityPlot.R")

argv <- commandArgs(T)
# argv[1] <- "/media/galaxy/My Passport/backup/olddata/03_ZHM/16s/shenglan/" 
# argv[2] <-  "20-25-30" 
# argv[3] <- "rare_test.pdf"
wdir<-argv[1]
rm_sample<-as.numeric(strsplit(argv[2],"-")[[1]])
fig_name <- argv[3]
##change 4870->5000
PL_func <- function(x) {
  ceiling(x/10^nchar(x)*10)*(10^(nchar(x)-1))
}

setwd(wdir)
outputFolder <- "4.figures/"
dat_raw<- read.table("1.data_for_statstics/otu_change_syringae.txt", row.names=1, sep="\t", header=T, blank.lines.skip = FALSE,stringsAsFactors = FALSE)
map<-read.table("1.data_for_statstics/map.tsv",row.names =1,sep="\t",stringsAsFactors = FALSE,header = T)
map <- map[-c(rm_sample),]
rownames(map) <- paste0("X",rownames(map))
Rhizo <- rownames(map[grepl("Rhizosphere",map$Description),])
root <- rownames(map[grepl("Root",map$Description),])
soil <- rownames(map[grepl("soil",map$Description),])
alldesign <- map
##
    dat_mutant <- dat_raw[,-dim(dat_raw)[2]]
    nsam <-dim(dat_mutant)[2] 
    dat_mutant <- dat_mutant[rowSums(dat_mutant) > 0, ]
    max_seq<- ceiling(sort(colSums(dat_mutant),decr=T)[1]/1000)*1000
    step_num <- max_seq/50
##make a file can be readable
    rownames(dat_mutant) <-gsub("^","COG", rownames(dat_mutant),perl=TRUE  )
    rownames(dat_mutant) <-gsub("^COGNew.ReferenceOTU","COG000", rownames(dat_mutant),perl=TRUE  )
    mid_table_name=paste("mutant_all",".txt",sep="")
    write.table(cbind("GeneFamily"=rownames(dat_mutant), dat_mutant), file=mid_table_name, sep="\t", row.names=FALSE, col.names=T, quote=FALSE)
    dat_mutant_GF <- readGeneFamilies(mid_table_name)
    dim(dat_mutant_GF$Data)
    rawplot <- emiel.diversityPlot.family(dat_mutant_GF, sample=1:nsam, max.sample.size=max_seq)
    aa <- matrix(0, 100,nsam )
    for(i in 1:nsam){aa[,i] <- rawplot[[i]]}
    lib_sizes <- colSums(dat_mutant)
    cols <- c("red","dimgrey","seagreen4","grey","steelblue2","tan1","seagreen","goldenrod","brown","mediumpurple3")
    names(cols)<-unique(map$Description)
    pdf(paste(outputFolder,fig_name,sep=""),width=7, height=4)
    max_y <- PL_func(round(max(aa[100,])))
    plot(c(0,max_seq), c(0,max_y), ylim=c(0,max_y), xlim=c(0,max_seq), frame.plot=F,col="transparent", ylab="OTUs detected", xlab="sequencing depth")
    par(mar=c(5,10,5,10),xpd=TRUE)
    colnames(aa) <-rownames(alldesign)
    for(i in Rhizo){
    ifelse(lib_sizes[i] < max_seq, x <- seq(1, lib_sizes[i], step_num), x <- seq(1, max_seq, step_num))
      points(aa[ 1:length(x),i], x=x, type="l",lty=4,lwd=0.7, col=cols[map[i,]$Description])
    }
    for(i in soil){
      ifelse(lib_sizes[i] < max_seq, x <- seq(1, lib_sizes[i], step_num), x <- seq(1, max_seq, step_num))
      points(aa[ 1:length(x),i], x=x, type="l",lty=1,lwd=0.7, col=cols[map[i,]$Description])
    }
    for(i in root){
      ifelse(lib_sizes[i] < max_seq, x <- seq(1, lib_sizes[i], step_num), x <- seq(1, max_seq, step_num))
      points(aa[ 1:length(x),i], x=x, type="l",lty=3,cex=0.4, lwd=0.7 ,col=cols[map[i,]$Description])
    } 
     legend(max_seq*0.8,max_y,legend = unique(map[Rhizo,]$Description),lty=4,col = cols[grep("Rhizo",names(cols))],bty="n",xpd = TRUE,cex=0.65,pt.cex = 1)
     legend(max_seq*0.8,max_y-1000,legend = unique(map[root,]$Description), lty=3,col = cols[grep("Root",names(cols))],bty="n",xpd = TRUE,cex=0.65,pt.cex = 1)
     legend(max_seq*0.8,max_y-2000,legend = unique(map[soil,]$Description), lty=1,col = cols[grep("Bulk",names(cols))],bty="n",xpd = TRUE,cex=0.65,pt.cex = 1)
     dev.off()


