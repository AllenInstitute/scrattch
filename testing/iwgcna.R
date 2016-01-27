###Code to run iterative WGCNA to identify cell clusters using high variance genes###
setwd("//aibsdata/rnaseqanalysis/SMARTer/SC_Core/analysis_1512/")

###Required libraries###
library(WGCNA)
library(limma)
library(flashClust)
library(matrixStats)
source("de.genes.R")
source("DESeq.var.R")
source("heatmap3.R")
source("cluster.R")
require(gplots)


###load data###
infoDir = "//aibsdata/rnaseqanalysis/SMARTer/LIMS_Alignment_Results//"
inDir = "//aibsdata/rnaseqanalysis/SMARTer/LIMS_Alignment_Results//"
outDir = "//aibsdata/rnaseqanalysis/SMARTer/SC_Core/analysis_1512/"
sample.info = file.path(infoDir, 'QCmetrics_RSC-005-9_LIMSupdate.csv')
ercc.file = file.path(inDir, 'ercc_counts_RSC-005-9_LIMS.csv')
rpkm.file = file.path(inDir, 'rpkm_RSC-005-9_LIMSgenesOnly.csv')
counts.file = file.path(inDir, 'counts_RSC-005-9_LIMSgenesOnly.csv')

plotsamp = 1;  #whether or not to plot variance curves from the Brennecke method
redogenes = 1;  #whether or not to recalculate high variance genes at each iteration
redopca = 0;   #whether or not to redo pca at each iteration
allrangevals = c(0,0.0625,0.25,1)   #variance thresholds to use to call high variance genes
mincountlimit = 200000; #minimum number of mapped reads (counts) to use


###load data###
##################
setwd(inDir)
sample = read.csv(sample.info, row.names=1, as.is=T)
ercc = read.csv(ercc.file, row.names=1, as.is=T)
ercc<-ercc[,-1]
rpkm = read.csv(rpkm.file, row.names=1, as.is=T)
counts = read.csv(counts.file, row.names=1,as.is=T)
ercc = ercc[,match(colnames(counts),colnames(ercc))]
rpkm = rpkm[,match(colnames(counts),colnames(rpkm))]
colnames(ercc) = gsub("\\.","-",colnames(ercc))
colnames(rpkm) = gsub("\\.","-",colnames(ercc))
colnames(counts) = gsub("\\.","-",colnames(ercc))
rownames(sample) = gsub("\\.","-",rownames(sample))
sample = sample[match(colnames(counts),rownames(sample)),]

setwd(outDir)
#create rda files to avoid loading large dataset each time
save(sample, file="sample.info.rda")
save(ercc, file="ercc.rda")
save(counts, file="counts.rda")
save(rpkm, file="rpkm.rda")

#load desired data files here
load(file="sample.info.rda")
load(file="ercc.rda")
load(file=paste(counts.file, "rda", sep="."))
load(file=paste(rpkm.file, "rda", sep="."))

####remove mitochondrial and tRNA genes###
genes = intersect(row.names(rpkm), row.names(counts))
select.genes = setdiff(genes, c(grep("rRNA",genes, value=T),grep("tRNA",genes, value=T),grep("^mt_",genes, value=T)))
all.cells = colnames(ercc)[colSums(ercc)>1000 & colSums(counts)>=mincountlimit]

counts1 = as.matrix(counts[select.genes,all.cells])
rpkm1 = as.matrix(rpkm[select.genes,all.cells])
ercc1 = as.matrix(ercc[, all.cells])
sample.all = sample[all.cells, ]
norm.dat.all = log2(rpkm1+1)

##set color scheme for plots
blue.red <-colorRampPalette(c("blue", "white", "red"))
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

wgcna_iteration<-function(subDir,clustnum=2,cells,min.cells,padj.th,prefix="top",minModuleSize,min.padj,cut_height) {
  dir.create(file.path(outDir, subDir), showWarnings=F)
  setwd(file.path(outDir, subDir))
  tmp=getGeneTree(norm.dat.all, ercc1, counts1, cells, min.cells, padj.th, prefix="top")
  geneTree = tmp[[1]]
  dissTOM= tmp[[2]]
  Vargenes = geneTree$labels
  write.csv(Vargenes, 'Vargenes.csv')
  
  tmp=cutGeneTree(norm.dat=norm.dat.all, select.cells=cells, geneTree=geneTree, dissTOM=dissTOM, minModuleSize=minModuleSize,min.cells=min.cells, padj.th=padj.th, min.padj=min.padj, cutHeight=cut_height)
  gene.mod= tmp[[1]]
  hc = tmp[[2]]
  gene.mod.val = tmp[[3]]
  
  markers = unlist(gene.mod)
  if ((is.null(gene.mod)) | (length(markers)<2)) {
    cl=rep(1,length(cells))
    names(cl)<-colnames(norm.dat.all[,cells])
    write.table(cl, file="cl.csv",quote=F)
    return(cl)
  }
  
  
  norm.dat = norm.dat.all[,cells]
  
  tmp.dat = norm.dat[markers,]
  print(length(markers))
  tmp.dat = tmp.dat - rowMeans(tmp.dat)
  tmp.dat1 = tmp.dat
  tmp.dat1[abs(tmp.dat1) < 0.5] = 0
  hc = hclust(dist(t(tmp.dat1)),method="ward")
  if (clustnum>0) {
    cl = cutree(hc, clustnum) 
  } else {
    tmp=cutDynamic(norm.dat,select.cells, hc, n.range=1:8,q=0.7) #final # of clusters will be determined within the range automatically based on the clustering results
    sc = sapply(tmp, function(x)sum(x$sc-100))
    tmp1 = tmp[[which.max(sc)]]
    cl = tmp1[[1]]
    de.df=tmp1[[2]]
    sc = tmp1[[3]]
  }
  
  cl = factor(cl)
  levels(cl) = 1:length(levels(cl))
  if (length(unique(cl))<2) {
    write.table(cl, file="cl.csv",quote=F)
    return(cl)
  }
  ##set side color bars here
  sample1 = sample.all[cells,]
  col1 = jet.colors(length(unique(cl)))[as.factor(cl)]
  col2 = jet.colors(length(unique(sample1$Batch)))[as.factor(sample1$Batch)]
  col3 = jet.colors(length(unique(sample1$Cre.line)))[as.factor(sample1$Cre.line)]
  col4 = jet.colors(length(unique(sample1$ROI)))[as.factor(sample1$ROI)]
  col5 = jet.colors(length(unique(sample1$transcriptomereads)))[as.factor(sample1$transcriptomereads)]
  #col4 = ifelse(sample1$tdTomato=='Y', 'red', 'gray')
  #col6 = jet.colors(length(unique(sample1$ERCC_percent)))[as.factor(sample1$ERCC_percent)]
  #col7 = jet.colors(length(unique(sample1$Collection_date)))[as.factor(sample1$Collection_date)]
  tmp.col = rbind(col1, col2, col3, col4, col5)
  markers=displayCl(norm.dat,cl, hc=hc, ColSideColors=tmp.col, prefix="tmp1",markers=markers)
  
  cl = cl[row.names(sample1)]
  cl.col = jet.colors(length(levels(cl)))[as.integer(cl)]
  colnames(tmp.col) = row.names(sample1)
  de.df = DE.genes.pw(norm.dat, paste("cl",cl, sep=""))
  markers=displayCl(norm.dat,cl, de.df=de.df, hc=hc,ColSideColors=tmp.col, prefix="top3",n.markers=3,labels=paste(cl, names(cl)),q=0.7)
  markers=displayCl(norm.dat,cl, de.df=de.df, hc=hc,ColSideColors=tmp.col, prefix="top10",n.markers=10,labels=paste(cl, names(cl)),q=0.7)
  for(x in names(de.df)){
    df = de.df[[x]]
    df = df[order(df$padj),]
    select = with(df, which(padj < 0.01 & abs(lfc) > 1))
    write.csv(df[select,],"de.genes.csv")
  }
  write.table(cl, file="cl.csv",quote=F)
  return(cl)
}


####first iteration###

clout<-wgcna_iteration(subDir='block0',clustnum=0,cells=all.cells,min.cells=5,padj.th=0.01,prefix="top",minModuleSize=5,min.padj=40,cut_height=0.995)

recursive_cluster(clout,'block0',keepcl=7:8);


recursive_cluster<-function(cl.set,headstring,keepcl=NULL) {
  allcls<-unique(cl.set)
  print(c(headstring,length(allcls)))
  if (length(allcls)>1) {
    if (!is.null(keepcl)) {
      allcls=keepcl;
    }
    for (cltyp in allcls) {
      if (length(which(cl.set==cltyp))>5) {
        headstring2=paste(headstring,cltyp,sep="_")
        print(c(headstring2,length(which(cl.set==cltyp))))
        clout<-wgcna_iteration(subDir=headstring2,clustnum=0,cells=names(cl.set)[cl.set==cltyp],min.cells=5,padj.th=0.01,prefix="top",minModuleSize=5,min.padj=40,cut_height=0.995)
        aa<-recursive_cluster(clout,headstring2)
      }
    }
  }
}



####gather clusters####
setwd(outDir)
clusterfolders<-c("0","0_1","0_1_1","0_2","0_2_3","0_2_4","0_2_4_2","0_2_4_2_3","0_2_6","0_2_8","0_3","0_4","0_5","0_5_1","0_5_3","0_6","0_6_5","0_7","0_8")
finalcl<-read.table(paste0("block",clusterfolders[1],"/cl.csv"),as.is=T,header=T)
for (ii in 2:length(clusterfolders)) {
  newcl<-read.table(paste0("block",clusterfolders[ii],"/cl.csv"),as.is=T,header=T)
  finalcl[rownames(newcl),1]<-paste(clusterfolders[ii],newcl[,1],sep="_")
}
write.csv(finalcl,file="iwgcna_clustering.csv")

####cluster comparison###
wgcnacl<-read.csv("iwgcna_clustering.csv",as.is=T)
seurcl<-read.csv("iseurat_clustering.csv",as.is=T)
seurcl[,1]<-gsub("\\.","-",seurcl[,1])
keepsampid<-suercl[,1]
for (ii in 1:nrow(seurcl)) {
  seurcl[ii,1]<-strsplit(seurcl[ii,1],"~")[[1]][1]
}
allcells<-matrix(0,nrow=ncol(rpkm),ncol=4);
rownames(allcells)<-colnames(rpkm)
colnames(allcells)<-c("w","s","int","intkeep")
allcells[wgcnacl[,1],1]<-wgcnacl[,2]
allcells[seurcl[,1],2]<-seurcl[,2]
crossclust<-table(allcells[,1],allcells[,2])
allcells[,3]<-paste(allcells[,1],allcells[,2],sep="v")
write.csv(crossclust,file="cross_cluster.csv")
ttt<-table(allcells[allcells[,3]!="0v0",3])
keepclusts<-names(ttt)[ttt>5]
allcells[,4]<-allcells[,3]
allcells[which(!(allcells[,4] %in% keepclusts)),4]<-"0"
write.csv(allcells,file="cluster_ids_w_and_s.csv")



