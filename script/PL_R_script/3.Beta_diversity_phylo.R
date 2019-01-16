###TIC rare 29000 sequences
#read Weighted UniFrac matrix of TIC
library(dendextend)
rm(list = ls())
argv <- commandArgs(T)
wdir<-argv[1]
input_wu <- argv[2]
outfig <- argv[3]


# argv[1]<-"/media/galaxy/My Passport/backup/olddata/03_ZHM/16s/shenglan/"
# argv[2]<-"3.data_for_figure/w_uni_tic/weighted_unifrac_group1_TIC_flt.txt"
# argv[3]<-"Beta_diversity_TIC_w_unifrac.pdf"
# argv[4]<-"5-6-7-8-9-20-25-30-35-36-37-38-39 "
# wdir <- "/media/galaxy/My Passport/backup/olddata/03_ZHM/16s/shenglan/"
# setwd(wdir)
system("mkdir 4.figures")
outfolder="./4.figures/"

rm_sample<-as.numeric(strsplit(argv[4],"-")[[1]])

W_unifrac<-read.table(input_wu,sep="\t",header = T,row.names = 1,stringsAsFactors = FALSE)
###set color
map<-read.table("1.data_for_statstics/map.tsv",sep = "\t",header = T,row.names = 1,stringsAsFactors = FALSE)
map <- map[-c(rm_sample),]
rownames(map) <- paste0("X",rownames(map))
mutant_name<-unique(map[!grepl("soil|Col",map$Description),]$Description)

soil<-rownames(map[grepl("soil",map$Description),])
soil_col<-rep("black",length(soil))
names(soil_col)<-soil

wt<-rownames(map[grepl("Col",map$Description),])
wt_col<-rep("blue",length(wt))
names(wt_col)<-wt
wt_col

mutant<-rownames(map[!grepl("soil|Col",map$Description),])
dup <- table(map[!grepl("soil|Col",map$Description),]$Description)
mutant_col<-rep(c("seagreen4","tan1","seagreen4","tan1"),dup)
names(mutant_col)<-mutant
mutant_col

W_unifrac_dist<-as.dist(W_unifrac)
W_unifrac_hclust<-hclust(W_unifrac_dist,method = "average")
col<-c(soil_col,wt_col,mutant_col)
col<-col[rownames(W_unifrac)[W_unifrac_hclust$order]]
col
rownames(W_unifrac)[W_unifrac_hclust$order]

###set pch
soil_pch<-rep(16,length(soil))
names(soil_pch)<-soil
map
rhizo<-rownames(map[grepl("Rhizosphere",map$Description),])
rhizo_pch<-rep(17,length(rhizo))
names(rhizo_pch)<-rhizo

root<-rownames(map[grepl("Root",map$Description),])
root_pch<-rep(18,length(root))
names(root_pch)<-root
pch<-c(soil_pch,rhizo_pch,root_pch)
pch<-pch[rownames(W_unifrac)[W_unifrac_hclust$order]]
pch


#pdf(paste(outfolder,"Beta_diversity_ACM_w_unifrac.pdf",sep=""), paper="special", width=10, height=5, horizontal = FALSE)
pdf(paste(outfolder,outfig,sep=""), width=16, height=5)
rownames(W_unifrac)[W_unifrac_hclust$order]
names(W_unifrac_hclust)
W_unifrac_hclust$labels<-""
W_unifrac_hclust$order
dend<-as.dendrogram(W_unifrac_hclust)
dend %>% set("leaves_pch", pch) %>%
  set("leaves_cex", 1.5) %>%
  set("leaves_col",col) %>%
  plot(main = "",ylab="Weight UNIFRAC distance",ylim=c(0,0.25),xlab="")
legend(30,0.25,legend =c("Soil","rhizosphere","root"),pch=c(16,17,18),bty ="n",xpd = TRUE,cex=0.9,pt.cex = 1.2)
legend(35,0.25,legend =c("wt", mutant_name),
       fill = c("blue","seagreen4","tan1","seagreen4","tan1")
       ,border = NA,bty ="n",xpd = TRUE,cex=0.9,pt.cex = 1.2,ncol = 1)
dev.off()

