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
load( "Rdata/cluster_viz_all_out.Rdata" )

#annotate clusters using gene accessibility scores
genes = read.delim("reference/gencode.vM16.gene.bed", sep="\t", header=F)

genes.gr = GRanges(genes[,1], 
    IRanges(as.numeric(genes[,2]), as.numeric(genes[,3])), name=genes[,4]
)

genes.tss.gr = 

marker.genes = c(
    "H2-Aa", "Cd74", "Cd80", "Aire",
    "Pigr", "Ly6a", "Vil1","Dclk1",
    "Cd3e", "Cd3d", "Cd3g", "Trac",
    "Ptprc", "Cd19", "Spib","Gfi1"
)

genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)]

# re-add the cell-by-bin matrix to the snap object
x.sp = addBmatToSnap(x.sp);
x.sp = createGmatFromMat(
    obj=x.sp, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=as.integer( args[4] )
)

# normalize the cell-by-gene matrix
x.sp = scaleCountMatrix(
    obj=x.sp, 
    cov=x.sp@metaData$passed_filters + 1,
    mat="gmat",
    method = "logRPM"
)

# smooth the cell-by-gene matrix using imputation
x.sp = runMagic(
    obj=x.sp,
    input.mat="gmat",
    step.size=3
)

pdf(paste( "figures/", args[2], "_", args[3], "_cluster_viz_marker_genes.pdf", sep="" ), height=13, width=10)

par(mfrow = c(4, 4))

for(i in 1:16){
    plotFeatureSingle(
        obj=x.sp,
        feature.value=x.sp@gmat[,i],
        method="umap", 
        main=marker.genes[i],
        point.size=0.1, 
        point.shape=19, 
        down.sample=10000,
        quantiles=c(0, 1)
		)
  }

dev.off()

save.image ( "Rdata/integrate_all_genes_out.Rdata" )