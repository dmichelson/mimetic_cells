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
load( "Rdata/filter_bmat_out.Rdata" )

#calculate diffusion maps for reducing dimensions
print( "snapatac pipeline: .. calculating diffusion maps ... ")
x.sp = runDiffusionMaps(
    obj=x.sp,
    input.mat="bmat", 
    num.eigs=50
)

#find number of significant components-- heuristic= stop when it looks like a blob
print( "snapatac pipeline: .. plotting pairwise dimensions ... ")
pdf(paste( "figures/", args[2], "_", args[3], "_dim_reduction_dimwise.pdf", sep="" ), height=4, width=5)
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
	
#initialize visualization - umap
print( "snapatac pipeline: .. initializing visualization, type=UMAP ... ")

x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:30, 
    method="umap",
    seed.use=10
)

#test clustering with different k
print( "snapatac pipeline: .. running knn/louvain clustering with variable k ... ")
pdf(paste( "figures/", args[2], "_", args[3], "_cluster_variable_k.pdf", sep="" ), height=10, width=10)
par(mfrow = c(3, 3))
for ( i in c(10,15,20,25,30,35,40,45,50) ) {
    print( paste0("checking k=",i) )

    x.sp = runKNN(
        obj=x.sp,
        eigs.dims=1:30,
        k=i
	)
	
    x.sp = runCluster(
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
        down.sample=20000,
        legend.add=FALSE
	    )
}
dev.off()

#graph based clustering using KNN, selecting the number of dims and k to use
print( "snapatac pipeline: .. running final clustering ... ")
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:30,
    k=30
    )
	
x.sp = runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10
    )

x.sp@metaData$cluster = x.sp@cluster

#viz data
print( "snapatac pipeline: .. visualizing data ... ")

x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:30, 
    method="umap",
    seed.use=10
    )

pdf(paste( "figures/", args[2], "_", args[3], "_viz_plain.pdf", sep="" ), height=5, width=6)
ggplot2::ggplot( as.data.frame(x.sp@umap), ggplot2::aes( x=x.sp@umap[,1], y=x.sp@umap[,2] ) ) +
    ggplot2::geom_point( color ="steelblue", alpha=0.5 ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP-1") +
    ggplot2::ylab("UMAP-2") +
    ggplot2::labs(col="gray")
dev.off()

pdf(paste( "figures/", args[2], "_", args[3], "_cluster_viz_qc.pdf", sep="" ), height=10, width=10)
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
print( "snapatac pipeline: .. visualizing each sample within data ... ")
pdf(paste( "figures/", args[2], "_", args[3], "_cluster_viz_per_sample.pdf", sep="" ), height=10, width=14)
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

print( "snapatac pipeline: .. saving work ... " )

save.image( "Rdata/cluster_viz_all_out.Rdata" )
