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
load( "Rdata/find_dars_out.Rdata" )

args = commandArgs(trailingOnly=TRUE)

#find motifs in each cluster's DARs using homer
print( "snapatac pipeline: .. running homer motif finding ... " )
system("which findMotifsGenome.pl")
for ( i in levels(x.sp@cluster) ) {
    motifs = runHomer(
	    x.sp[,idy.ls[[as.integer(i)]],"pmat"], 
    	mat = "pmat",
    	path.to.homer = "/n/scratch2/dam41/homer/bin/findMotifsGenome.pl",
    	result.dir = paste0("peaks/homer/",args[3],"/cluster_",i),
    	num.cores=args[4],
    	genome = 'mm10',
    	motif.length = 10,
    	scan.size = 300,
    	optimize.count = 2,
    	background = 'automatic',
    	local.background = FALSE,
    	only.known = TRUE,
    	only.denovo = FALSE,
    	fdr.num = 5,
    	cache = 100,
    	overwrite = TRUE,
    	keep.minimal = FALSE,
    )
}

print( "snapatac pipeline: .. saving work ... " )

save.image( "Rdata/run_homer_motif.Rdata" )