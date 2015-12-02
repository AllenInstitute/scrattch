library(dplyr)
# ercc_counts <- read.csv("//AIBSData2/mct-t200/HiSeq/Analysis/032315_hiseq/1679_cells_passQC/ercc_counts.csv",row.names=1)
# ercc_counts <- ercc_counts[!grepl("tdT",rownames(ercc_counts)),]
# ercc_count_mat <- as.matrix(ercc_counts,dimnames=c(rownames(ercc_counts),names(ercc_counts)))
# 
# gene_counts <- read.csv("//AIBSData2/mct-t200/HiSeq/Analysis/032315_hiseq/1679_cells_passQC/counts.csv",row.names=1)
# gene_counts <- gene_counts[!grepl("^mt|^c",rownames(gene_counts)),]
# gene_count_mat <- as.matrix(gene_counts,dimnames=c(rownames(gene_counts),names(gene_counts)))
# 
# gene_rpkm <- read.csv("//AIBSData2/mct-t200/HiSeq/Analysis/032315_hiseq/1679_cells_passQC/rpkm.csv",row.names=1)
# gene_rpkm <- gene_rpkm[!grepl("^mt|^c",rownames(gene_rpkm)),]
# gene_rpkm_mat <- as.matrix(gene_rpkm,dimnames=c(rownames(gene_rpkm),names(gene_rpkm)))

run_brennecke <- function(ercc_count_mat, 
                          gene_count_mat, 
                          minBiolDisp = 1, 
                          plotsamp = FALSE) {
  library(DESeq)
  library(statmod)
  library(matrixStats)
  
  # Remove tdTomato from the ERCC table
  ercc_count_mat <- ercc_count_mat[!grepl("tdTomato",rownames(ercc_count_mat)),]
  
  # Estimate size factors for ERCCs and genes
  sf_ercc <- estimateSizeFactorsForMatrix(ercc_count_mat)
  sf_gene <- estimateSizeFactorsForMatrix(gene_count_mat)
  
  # Adjust counts by dividing by the size factor
  nCounts_ercc <- t( t(ercc_count_mat) / sf_ercc)
  nCounts_gene <- t( t(gene_count_mat) / sf_gene)
  
  # Calculate means, variances, and cv^2 for each ercc
  means_ercc <- rowMeans( nCounts_ercc )
  vars_ercc <- rowVars(nCounts_ercc)
  cv2_ercc <- vars_ercc / (means_ercc**2)
  
  # Calculate the minimum mean for the fit as the 80th percentile of means for which cv^2 is > 0.3
  minMeanForFit <- unname(quantile(means_ercc[which(cv2_ercc > 0.3)], 0.8))
  
  # Determine which ERCCs have means above the minMeanForFit, and perform 
  # generalized linear model fitting using those ERCCs.
  useForFit <- means_ercc >= minMeanForFit
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means_ercc[useForFit] ),cv2_ercc[useForFit] )
  
  # Extract variables from the fit results for later use.
  xi <- mean( 1 / sf_ercc )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"] - xi )
  print(paste0("a0 = ",round(a0,5),"; a1 = ",round(a1,5)))
    
  ##### test genes against the ERCC fit ####
  means_gene <- rowMeans( nCounts_gene )
  vars_gene <- rowVars(nCounts_gene)
  cv2_gene <- vars_gene / means_gene^2
  psia1theta <- mean( 1 / sf_gene ) + a1 * mean( sf_ercc / sf_gene )
  m <- nrow(gene_count_mat)
  cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
  testDenom <- ( means_gene * psia1theta + means_gene^2 * cv2th ) / ( 1 + cv2th/m )
  p <- 1 - pchisq( vars_gene * (m-1) / testDenom, m-1 )
  padj <- p.adjust( p, "BH" )
  sig <- padj < .05
  sig[is.na(sig)] <- FALSE
  print(paste0("Significant Genes: ",sum(sig)))
  
  #####plot samples####
  #windows();
  if (plotsamp) {
    df <- ncol(ercc_count_mat) - 1
    plot(NULL, xaxt="n", yaxt="n",log="xy", xlim = c( 1e-3, 1e4 ), cex.axis=2,cex.lab=2,ylim = c( .005, ncol(gene_count_mat) ),xlab = "average normalized UMI count", ylab = "squared coefficient of variation (CV^2)" )
    axis(1, 10^(-3:4), c( "0.001","0.01","0.1", "1", "10", "100", "1000",expression(10^4)),cex.axis=1.5 )
    axis(2, 10^(-2:2), c( "0.01", "0.1", "1", "10","100" ), las=2,cex.axis=1.5 )
    abline( h=10^(-2:2), v=10^(-1:5), col="#D0D0D0", lwd=2 )
    # Plot the plant genes, use a different color if they are highly variable
    points( means_gene, cv2_gene, pch=20, cex=.2,col = ifelse( padj < .05, "#C0007090", "darkgreen" ) )
    points( means_ercc, cv2_ercc, pch=20,cex=1, col="blue")
    #Add the technical noise fit, as before
    xg <- 10^seq( -2, 6, length.out=1000 )
    lines( xg, (xi+a1)/xg + a0, col="#FF000080", lwd=3 )
    # Add a curve showing the expectation for the chosen biological CV^2 thershold
    lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .975, df ) / df , col="grey", lwd=2, lty="dashed" )
    lines( xg, ( (xi+a1)/xg + a0 ) * qchisq( .025, df ) / df , col="grey", lwd=2, lty="dashed" )
    lines( xg, psia1theta/xg + a0 + minBiolDisp, lty="dashed", col="#C0007090", lwd=3 )
  }
  
  ####Table of genes####
  log2RelExpr_gene <- log2(nCounts_gene / means_gene)
  colnames(log2RelExpr_gene) <- 1:ncol(log2RelExpr_gene)
  highVarTable <- data.frame(geneSymbol = rownames(gene_count_mat)[sig], 
                             meanNormCount = means_gene[sig],
                             strongest = factor(colnames(log2RelExpr_gene)[apply(log2RelExpr_gene[sig,], 1, which.max)]),
                             log2RelExpr_gene[sig,],
                             check.names=FALSE)
  fullVarTable <- data.frame(geneSymbol = rownames(gene_count_mat),
                             pval = padj,
                             meanNormCount = means_gene,
                             cv2 = cv2_gene, 
                             log2RelExpr_gene,
                             check.names=FALSE)
  outlist <- c(ref = list(nCounts_ercc),
               samp = list(nCounts_gene),
               hi = list(highVarTable),
               full = list(fullVarTable))
  return(outlist)
}

block_cluster <- function(titlename, 
                          ercc_count_mat, gene_count_mat, minBiolDisp, 
                          gene_rpkm_mat, samples, popcolvec, keeplim, redopca = FALSE, textvec, meth, mindist) {

  brennecke_results <- run_brennecke(ercc_count_mat,gene_count_mat,minBiolDisp)
  
  sig_genes <- rownames(brennecke_results$hi)
  sig_gene_rpkm_mat <- gene_rpkm_mat[rownames(gene_rpkm_mat) %in% sig_genes, samples]
  
  if (redopca) {
    
    collapse_results <- collapse_genes(sig_gene_rpkm_mat[rowSums(sig_gene_rpkm_mat) > 0,], 5, 0)
    sig_gene_rpkm_mat <- collapse_results[['collapsed']]
    colnames(sig_gene_rpkm_mat) <- colnames(sig_gene_rpkm_mat)
    pca1 <- prcomp(scale(t(log2(sig_gene_rpkm_mat[collapse_results$clustcoh_meanr < 1,] + 1))))

  } else {
    
    sig_gene_rpkm_mat <- sig_gene_rpkm_mat[which(apply(sig_gene_rpkm_mat, 1, max) >= 0),]
    pca1 <- prcomp(scale(t(log2(sig_gene_rpkm_mat + 1))))
  
  }
  
  plot(pca1$x[,1:2],pch=16,col=popcolvec[1,colnames(ndat) %in% block1],main=titlename)
  plot(pca1$x[,1:2],pch='',col=popcolvec[1,colnames(ndat) %in% block1],main=titlename)
  text(pca1$x[,1:2],textvec[colnames(ndat) %in% block1],col=popcolvec[1,colnames(ndat) %in% block1])
  if (keeplim==0) {
    keepers<-sigcomps(scale(t(log2(startmat[outlist$clustcoh_meanr<1,]+1))))
  } else {
    keepers<-keeplim;
  }
  print(head(summary(pca1)$importance[2,]));
  print(keepers);
  keepgenes<-c();
  newdistmat<-matrix(0,nrow=nrow(pca1$x),ncol=nrow(pca1$x));
  coordmat<-matrix(0,nrow=nrow(pca1$x),ncol=keepers)
  for (ii in 1:keepers) {
    tempdist<-outer(pca1$x[,ii],pca1$x[,ii],"-");
    if (meth=='ward') {
      newdistmat<-newdistmat+(abs(tempdist)*(summary(pca1)$importance[2,ii]))^2;  
    } else {
      newdistmat<-newdistmat+abs(tempdist)*(summary(pca1)$importance[2,ii]);
    }
    coordmat[,ii]<-pca1$x[,ii]*summary(pca1)$importance[2,ii];
    if (nrow(pca1$rotation)>50) {
      keepgenes<-c(keepgenes,rownames(pca1$rotation)[order(-abs(pca1$rotation[,ii]))[1:50]])
    } else {
      keepgenes<-c(keepgenes,rownames(pca1$rotation)[order(-abs(pca1$rotation[,ii]))])
    }
  }
  keepgenes<-unique(keepgenes)
  rc3=as.dendrogram(hclust(as.dist(newdistmat),method=meth));
  rc4=rc3;
  temp<-cmdscale(as.dist(newdistmat),k=3,eig=T);
  rownames(temp$points)<-rownames(newdistmat);
  titlename<-paste(titlename, ", with ",keepers," PCs and ",nrow(startmat)," genes",sep='')
  if (mindist!=0) {
    breakvals<-((0:100)*(mindist/100))-mindist;
    heatmap.3(as.matrix(-newdistmat),main=titlename,Rowv=rc4,Colv=rc4,breaks=breakvals,RowSideColors=popcolvec[2,colnames(ndat) %in% block1],ColSideColors=popcolvec[1,colnames(ndat) %in% block1],trace='none',scale='none',col=jet.colors,hclustfun=hclust2,distfun=dist2p);  
  } else {
    heatmap.3(as.matrix(-newdistmat),main=titlename,Rowv=rc4,Colv=rc4,RowSideColors=popcolvec[2,colnames(ndat) %in% block1],ColSideColors=popcolvec[1,colnames(ndat) %in% block1],trace='none',scale='none',col=jet.colors,hclustfun=hclust2,distfun=dist2p);  
  }
  breakvals2<-((0:400)*0.01)-2;
  heatmap.2(as.matrix(log2(startmat[rownames(startmat) %in% keepgenes,]+1)),trace='none',Colv=rc4,hclustfun=hclust2,distfun=dist2p,col=bluered,breaks=breakvals2,scale='row',ColSideColors=popcolvec[1,colnames(ndat) %in% block1],cexRow=0.2,cexCol=0.2)
  temprcblock1=hclust(as.dist(newdistmat),method=meth)
  outlist2<-list()
  outlist2$rcout=temprcblock1;
  outlist2$dend=rc4;
  outlist2$mds<-temp;
  outlist2$distmat<-newdistmat;
  outlist2$pcs<-summary(pca1)$importance[2,];
  outlist2$pca1<-pca1;
  outlist2$coord<-coordmat;
  return(outlist2);
}