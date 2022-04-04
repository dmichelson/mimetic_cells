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
load( "Rdata/remove_clusters_out.Rdata" )

#add cell-by-bin matrix at desired resolution
print( "snapatac pipeline: .. re-adding bmat ... ")
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=as.integer( args[4] ))

#binarize the cell-by-bin matrix
print( "snapatac pipeline: .. binarizing bmat ... ")
x.sp = makeBinary(x.sp, mat="bmat")

#exclude blacklisted regions
print( "snapatac pipeline: .. excluding blacklist ... ")
black_list = read.table("reference/mm10.blacklist.bed.gz")
black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
	)
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr))

if(length(idy) > 0)
	{
	x.sp = x.sp[,-idy, mat="bmat"]
	}
x.sp

#exclude unwanted chromosomes, i.e random, ChrM
print( "snapatac pipeline: .. excluding unwanted chromosomes ... ")
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM|chrUn|chrY", seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0)
	{
	x.sp = x.sp[,-idy, mat="bmat"]
	}
x.sp

#exclude top 5% bins that overlap with invariant features like housekeeping genes
print( "snapatac pipeline: .. excluding top 5% bins ... ")
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]
x.sp

#filter low bin coverage cells (<1000) -- optional but recommended
print( "snapatac pipeline: .. excluding low coverage cells ... ")
idx = which(Matrix::rowSums(x.sp@bmat) > 1000);
x.sp = x.sp[idx,];
x.sp

#rerun clustering
print( "snapatac pipeline: .. running diffusion maps ... ")
x.sp = runDiffusionMaps(
    obj=x.sp,
    input.mat="bmat", 
    num.eigs=50
)

#revisualize the data
print( "snapatac pipeline: .. initiating viz ... ")
x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:20, 
    method="umap",
    seed.use=10
)

#refind number of significant components-- heuristic= stop when it looks like a blob
print( "snapatac pipeline: .. plotting dimwise reduction ... ")
pdf(paste( "figures/", args[2], "_", args[3], "_dim_reduction_dimwise_postCluster.pdf", sep="" ), height=4, width=5)
plotDimReductPW(
    obj=x.sp, 
    eigs.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
    )
dev.off()
	
#rerun clustering with different k
print( "snapatac pipeline: .. running knn/louvain clustering with variable k ... ")
pdf(paste( "figures/", args[2], "_", args[3], "_cluster_variable_k_postCluster.pdf", sep="" ), height=10, width=10)
print( paste0("snapatac pipeline: .. k=",i," ... ") )
par(mfrow = c(3, 3))
for ( i in c(10,15,20,25,30,35,40,45,50) ) {
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:20,
    k=i
	)
	
x.sp=runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10
	)

x.sp@metaData$cluster = x.sp@cluster
plotViz(
    obj=x.sp,
    method="umap", 
    main=paste("MEC Cluster, k=",i,sep=""),
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
    down.sample=5000,
    legend.add=FALSE
	)
}
dev.off()

#graph based clustering using KNN, selecting the number of dims and k to use
print( "snapatac pipeline: .. running clustering with desired k ... ")

x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:20,
    k=45
	)
	
x.sp = runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10
	)

x.sp@metaData$cluster = x.sp@cluster

x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:20, 
    method="umap",
    seed.use=5
)

pdf(paste( "figures/", args[2], "_", args[3], "_viz_postRemoval_plain.pdf", sep="" ), height=5, width=6)
ggplot2::ggplot( as.data.frame(x.sp@umap), ggplot2::aes( x=x.sp@umap[,1], y=x.sp@umap[,2] ) ) +
    ggplot2::geom_point( color ="steelblue", alpha=0.5 ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP-1") +
    ggplot2::ylab("UMAP-2") +
    ggplot2::labs(col="gray")
dev.off()

print( "snapatac pipeline: .. visualizing clusters ... ")

pdf(paste( "figures/", args[2], "_", args[3], "_cluster_viz_postCluster_qc.pdf", sep="" ), height=10, width=10)
par(mfrow = c(2, 2))
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
    down.sample=20000,
    legend.add=FALSE
)
	
plotFeatureSingle(
    obj=x.sp,
    feature.value=log(x.sp@metaData[,"passed_filters"]+1,10),
    method="umap", 
    main="Merged MEC Read Depth",
    point.size=0.2, 
    point.shape=19, 
    down.sample=20000,
    quantiles=c(0.01, 0.99)
)
	
plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@metaData$peak_region_fragments / x.sp@metaData$passed_filters,
    method="umap", 
    main="Merged MEC FRiP",
    point.size=0.2, 
    point.shape=19, 
    down.sample=20000,
    quantiles=c(0.01, 0.99) # remove outliers
)
	
plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@metaData$duplicate / x.sp@metaData$total,
    method="umap", 
    main="Merged MEC Duplicate",
    point.size=0.2, 
    point.shape=19, 
    down.sample=20000,
    quantiles=c(0.01, 0.99) # remove outliers
)

dev.off()

#look at distribution of individual samples within larger clustering
print( "snapatac pipeline: .. visualizing samples individually ... ")

pdf(paste( "figures/", args[2], "_", args[3], "_cluster_viz_postCluster_per_sample.pdf", sep="" ), height=10, width=14)
par(mfrow = c(2, 3))
plotViz(
    obj= x.sp,
    method="umap", 
    main="Merged MEC Sample",
    point.size=0.2, 
    point.shape=19, 
    point.color=x.sp@sample, 
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=20000,
    legend.add=TRUE
)

plotViz(
    obj= x.sp[order(x.sp@sample=="WT1")],
    method="umap", 
    main="WT1 MEC Sample",
    point.color=x.sp[order(x.sp@sample=="WT1")]@sample=="WT1", 
    point.size=0.2, 
    point.shape=19,
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=20000,
    legend.add=TRUE
)

plotViz(
    obj= x.sp[order(x.sp@sample=="WT2")],
    method="umap", 
    main="WT2 MEC Sample",
    point.color=x.sp[order(x.sp@sample=="WT2")]@sample=="WT2", 
    point.size=0.2, 
    point.shape=19,
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=20000,
    legend.add=TRUE
)

plotViz(
    obj= x.sp[order(x.sp@sample=="KO1")],
    method="umap", 
    main="KO1 MEC Sample",
    point.color=x.sp[order(x.sp@sample=="KO1")]@sample=="KO1", 
    point.size=0.2, 
    point.shape=19,
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=20000,
    legend.add=TRUE
)

plotViz(
    obj= x.sp[order(x.sp@sample=="KO2")],
    method="umap", 
    main="KO2 MEC Sample",
    point.color=x.sp[order(x.sp@sample=="KO2")]@sample=="KO2", 
    point.size=0.2, 
    point.shape=19,
    text.add=FALSE,
    text.size=1.5,
    text.color="black",
    down.sample=20000,
    legend.add=TRUE
)

dev.off()

print( "snapatac pipeline: .. saving work ... ")

save.image( "Rdata/cluster_viz_subset_out.Rdata" )