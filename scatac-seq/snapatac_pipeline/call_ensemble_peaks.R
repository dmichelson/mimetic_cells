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
load( "Rdata/cluster_clusters_out.Rdata" )

#call ensemble peaks
system("which snaptools")
system("which macs2")

#make ensemble peaks for all clusters
print( "snapatac pipeline: .. calling ensemble peaks with MACS2 ... ")

# call peaks for all cluster with more than 100 cells
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 100)]
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i])
    runMACS(
        obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
        output.prefix=paste0("peaks/macs2/",args[2],"_",args[3],"_peaks", gsub(" ", "_", clusters.sel)[i]),
        path.to.snaptools="/home/dam41/.conda/envs/snaptools_env/bin/snaptools",
        path.to.macs="/n/app/macs2/2.1.1.20160309/bin/macs2",
        gsize="mm", # mm, hs, etc
        buffer.size=500, 
        num.cores=as.integer(args[4]),
        macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
        tmp.folder=tempdir()
   )
 }, mc.cores=as.integer(args[4])
)
 
# assuming all .narrowPeak files in the current folder are generated from the clusters
print( "snapatac pipeline: .. finding peak files ... ")
peaks.names = system("ls peaks/macs2 | grep narrowPeak", intern=TRUE)
peak.gr.ls = lapply(paste0("peaks/macs2/",peaks.names), function(x){
    peak.df = read.table(x)
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls))
peak.gr

#save combined peaks file
print( "snapatac pipeline: .. merging peak files ... ")
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = paste0("peaks/combined_peaks/",args[2],"_",args[3],"_peaks_combined.bed"),append=FALSE,
	quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
	row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
	fileEncoding = ""
)

print( "snapatac pipeline: .. saving work ... ")

save.image( "Rdata/call_ensemble_peaks_out.Rdata" )