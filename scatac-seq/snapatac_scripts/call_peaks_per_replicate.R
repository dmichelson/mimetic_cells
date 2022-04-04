#!/usr/bin/env Rscript

#snapatac Rscript
#feed inputs as "Rscript run_snapatac_generic.R *wd *snap_file_name *date *n_cores [ *arg6...]"

#outside the environment, load homer and macs2 into path. Homer is in scratch directory and may need to be reinstalled...

#set up environment
args = commandArgs(trailingOnly=TRUE)

.libPaths(.libPaths()[2])

setwd(args[1])
getwd()

library("tidyverse")
library("SnapATAC")
library("viridisLite")
library("GenomicRanges")
library("chromVAR")
library("motifmatchr")
library("SummarizedExperiment")
library("BSgenome.Mmusculus.UCSC.mm10")
library("JASPAR2016")
library("BiocParallel")
library("pheatmap")
library("cluster")
library("BuenColors")

load( "Rdata/subclustered_postaire_out.Rdata" )

args = c( getwd(), "WT1_WT2_KO1_KO2_merged", "2020-11-25", "8" )

custom_colormaps = jdb_color_maps
names(custom_colormaps) = seq(1:17)

#run macs2
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 0)]
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i])
    runMACS(
        obj=x.sp[which(x.sp@cluster==clusters.sel[i] & x.sp@sample %in% c("WT1","KO1")),], 
        output.prefix=paste0("peaks/",args[3],"/",args[2],"_",args[3],"_peaks", gsub(" ", "_", clusters.sel)[i],"_","rep1"),
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

clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 0)]
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i])
    runMACS(
        obj=x.sp[which(x.sp@cluster==clusters.sel[i] & x.sp@sample %in% c("WT2","KO2")),], 
        output.prefix=paste0("peaks/",args[3],"/",args[2],"_",args[3],"_peaks", gsub(" ", "_", clusters.sel)[i],"_","rep2"),
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

