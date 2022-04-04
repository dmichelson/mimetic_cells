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
load( "Rdata/preprocess_snap_out.Rdata" )

#add cell-by-bin matrix at desired resolution
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=as.integer( args[4] ))

#binarize the cell-by-bin matrix
x.sp = makeBinary(x.sp, mat="bmat")

#exclude blacklisted regions
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
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM|chrUn|chrY", seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0)
	{
	x.sp = x.sp[,-idy, mat="bmat"]
	}
x.sp

#exclude top 5% bins that overlap with invariant features like housekeeping genes
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)

pdf(paste( "figures/", args[2], "_", args[3], "_bin_coverage_hist.pdf", sep="" ), height=4, width=5)
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(Bin Cov)", 
    col="steelblue", 
    xlim=c(0, 5)
	)
abline( v=bin.cutoff, col="red", lty=2 )
abline( v=1, col="blue", lty=2 )
dev.off()
	
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]
x.sp

#filter low coverage cells (<1000) -- optional but recommended
idx = which(Matrix::rowSums(x.sp@bmat) > 1000);
x.sp = x.sp[idx,];
x.sp

print( "snapatac pipeline: .. saving work ... " )

save.image( "Rdata/filter_bmat_out.Rdata")