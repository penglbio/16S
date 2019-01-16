rm(list = ls())
library("Cairo")
argv <- commandArgs(T)
argv[1] <- "/media/galaxy/My Passport/backup/olddata/03_ZHM/16s/shenglan/" 
argv[2] <-  "5-6-7-8-9-20-25-30-35-36-37-38-39" 
  

wdir<-argv[1]
rm_sample<-as.numeric(strsplit(argv[2],"-")[[1]])


##map file
setwd(wdir)
map<-read.table("3.data_for_figure/map.tsv",header=T,row.names=1,sep = "\t")
map <- map[-c(rm_sample),]
rownames(map) <- paste0("X",rownames(map))
Rhizo <- rownames(map[grepl("Rhizosphere",map$Description),])
root <- rownames(map[grepl("Root",map$Description),])
soil <- rownames(map[grepl("soil",map$Description),])

# ----------------------------------------------------------phylaRA---------------------------------------------------------
##phylum_RA soil rhizo root
phyla_RA <- read.table("../1.data_for_stastic/Phy_RA_table.txt", row.names=1, sep="\t", header=T, blank.lines.skip = FALSE,stringsAsFactors = FALSE)
bs_phyla<-phyla_RA[,soil]
rhizo_phyla<-phyla_RA[,Rhizo]
root_phyla<-phyla_RA[,root]
col_fuls<-c("aquamarine","aquamarine4","blue","blue4","blueviolet","brown","brown1","brown4","darkgreen",
            "darkolivegreen1","darkorange","darkorange2","darkslateblue","deeppink","deeppink3","gray48",
            "lightsalmon","palevioletred1","yellow","azure1","azure4")
##constrite_categrige
barstack_cat <- function(con_cat,phyla_RA){
  # con_cat <- root_phyla
  con_cat_names<-rownames(phyla_RA[with(phyla_RA, order(-phyla_RA[,ncol(phyla_RA)])), ])
  cat_names <-as.vector(na.exclude(con_cat_names[!grepl("not_assigned",con_cat_names)][1:20]))
  cat_cols <- col_fuls[1:length(cat_names)]
  names(cat_cols) <- cat_names
  other_cat<-setdiff(con_cat_names,cat_names)
  cols_other<-rep("black",length(other_cat))
  names(cols_other)<-other_cat
  con_catcols_legend<-c(cat_cols,cols_other)
  con_catcols_legend[rownames(con_cat)]
  con_cat<-con_cat[con_cat_names,]
  list(con_catcols_legend,con_cat)
}

stackbarplot <- function(vpar,con,classify,concat,xlend1,cpch,main){
  par(mar=vpar,xpd=TRUE)
  con_family <- concat[[2]]
  colnames(con_family)<-map[con,]$Description
  bp<-barplot(as.matrix(con_family), main="", lwd=1, col=concat[[1]], cex.names=.7, ylab = "RA [\u2030]", las=2, axisnames =F,border = NA)
  reps <- table(colnames(con_family))
  reps <- reps[unique(colnames(con_family))]
  nsample<- length(reps)
  points(bp,rep(-50,ncol(con_family)),pch=cpch,col=rep(col_fuls[1:nsample],reps))
  text(8,1100,paste("Taxonomy of",main,"ACM /", classify,sep = " "),font = 2,cex=1.1,xpd=TRUE)
  legend(xlend1,-100,legend=names(reps),pch=cpch,col=col_fuls[1:nsample],bty="n",cex=0.6,pt.cex = 0.5,ncol=2)
}



bs_catname_phy<- barstack_cat(bs_phyla,phyla_RA)
rhizo_catname_phy<- barstack_cat(rhizo_phyla,phyla_RA)
root_catname_phy<- barstack_cat(root_phyla,phyla_RA)

##plot phylum
CairoPDF("figure6_taxonomy/Phyla_Taxonomy_ACM_new.pdf",width = 10,height = 6)
layout(matrix(c(1,2,3,4),2,2),widths=c(3,1))
stackbarplot(vpar=c(4,4,4,8),Rhizo,"Phyrum",rhizo_catname_phy,4,17,"rhizo")
stackbarplot(vpar=c(5,4,4,8),root,"Phyrum",root_catname_phy,4,16,"root")
stackbarplot(vpar=c(3,4,4,8.5),soil,"Phyrum",bs_catname_phy,1,15,"soil")
legend(17,1000,legend = c(names(root_catname_phy[[1]][!grepl("black",root_catname_phy[[1]])]),"other classify") ,
       col =unique(root_catname_phy[[1]]), pch=16,bty = "n",xpd = TRUE,cex=0.5,pt.cex = 0.8)
dev.off()
# -------------------------------------------------------Family------------------------------------------------------------------------------------
family_RA <- read.table("3/Family_RA_table.txt", row.names=1, sep="\t", header=T, blank.lines.skip = FALSE,stringsAsFactors = FALSE)
bs_family<-family_RA[,soil]
rhizo_family<-family_RA[,Rhizo]
root_family<-family_RA[,root]

bs_catname_fam<- barstack_cat(bs_family,family_RA)
rhizo_catname_fam<- barstack_cat(rhizo_family,family_RA)
root_catname_fam<- barstack_cat(root_family,family_RA)

##plot family
CairoPDF("4.figures/family_Taxonomy_ACM.pdf",width = 10,height = 6)
layout(matrix(c(1,2,3,4),2,2),widths=c(2.6,1))
stackbarplot(vpar=c(4,4,4,8),Rhizo,"family",rhizo_catname_fam,4,17,"rhizo")
stackbarplot(vpar=c(5,4,4,8),root,"family",root_catname_fam,4,16,"root")
stackbarplot(vpar=c(3,4,4,8.5),soil,"family",bs_catname_fam,1,15,"soil")
legend(17,1000,legend = c(names(root_catname_fam[[1]][!grepl("black",root_catname_fam[[1]])]),"other classify") ,
       col =unique(root_catname_fam[[1]]), pch=16,bty = "n",xpd = TRUE,cex=0.5,pt.cex = 0.8)
dev.off()
