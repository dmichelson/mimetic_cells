#!/usr/bin/env Rscript

#snapatac Rscript
#feed inputs as "Rscript run_snapatac_generic.R *wd *snap_file_name *date *n_cores [ *arg6...]"

#outside the environment, load homer and macs2 into path. Homer is in scratch directory and may need to be reinstalled...

#set up environment
args = commandArgs(trailingOnly=TRUE)

.libPaths(.libPaths()[2])

setwd(args[1])
print( "snapatac pipeline: .. the working directory is ... ")
getwd()

suppressMessages( suppressWarnings( library("tidyverse") ) )
suppressMessages( suppressWarnings( library("SnapATAC") ) )
suppressMessages( suppressWarnings( library("viridisLite") ) )
suppressMessages( suppressWarnings( library("GenomicRanges") ) )
suppressMessages( suppressWarnings( library("chromVAR") ) )
suppressMessages( suppressWarnings( library("motifmatchr") ) )
suppressMessages( suppressWarnings( library("SummarizedExperiment") ) )
suppressMessages( suppressWarnings( library("BSgenome.Mmusculus.UCSC.mm10") ) )
suppressMessages( suppressWarnings( library("JASPAR2016") ) )
suppressMessages( suppressWarnings( library("BiocParallel") ) )

print( "snapatac pipeline: .. loading data ... ")
load( "Rdata/refilter_bmat_out.Rdata" )

# hierarchical clustering
# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
	SnapATAC::colMeans(x.sp[x,], mat="bmat")
} )

# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");

pdf(paste( "figures/", args[2],"_",args[3], "_hclust_viz.pdf", sep="" ), height=10, width=10)
plotViz(
    obj=x.sp,
    method="umap", 
    main="Merged MEC Cluster",
    point.color=x.sp@cluster, 
    point.size=1, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=TRUE,
    text.size=1.5,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    legend.add=FALSE
	)

plot(hc, hang=-1, xlab="")
dev.off()

save.image( "Rdata/cluster_clusters_out.Rdata" )