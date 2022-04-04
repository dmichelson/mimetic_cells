#!/usr/bin/env Rscript

#snapatac Rscript
#feed inputs as "Rscript run_snapatac_generic.R *wd *snap_file_name *date *n_cores [ *arg6...]"

#outside the environment, load homer and macs2 into path. Homer is in scratch directory and may need to be reinstalled...

#set up environment
args = commandArgs(trailingOnly=TRUE)
args = c( getwd(), "WT1_WT2_KO1_KO2_merged", "2021-04-03", "8" )

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

x.sp <- readRDS("Rdata/subclustered_postaire_out_2021-02-16.rds")

for (i in levels( unique( x.sp@cluster ) ) ) {
    print( paste0("Total cells in cluster ",i, " = ", sum( ( x.sp@cluster == i ) ) ) )
    print( paste0("WT cells in cluster ",i, " = ", sum( ( x.sp@sample %in% c( "WT1", "WT2" ) ) & ( x.sp@cluster == i ) ) ) )
    print( paste0("KO cells in cluster ",i, " = ", sum( ( x.sp@sample %in% c( "KO1", "KO2" ) ) & ( x.sp@cluster == i ) ) ) )
    print( paste0("WT fraction in cluster ",i, " = ", sum( ( x.sp@sample %in% c( "WT1", "WT2" ) ) & ( x.sp@cluster == i ) ) / sum( x.sp@cluster == i ) ) )
    print( paste0("KO fraction in cluster ",i, " = ", sum( ( x.sp@sample %in% c( "KO1", "KO2" ) ) & ( x.sp@cluster == i ) ) / sum( x.sp@cluster == i ) ) )
}

write_delim( as.data.frame(x.sp@sample), "analysis/2021-04-03_cluster_composition/x_sp_sample_2021-04-03.txt", delim="\t", col_names=F )
write_delim( as.data.frame(x.sp@cluster), "analysis/2021-04-03_cluster_composition/x_sp_cluster_2021-04-03.txt", delim="\t", col_names=F )

#calculate sample numbers per cluster

cluster_matrix = matrix( 0, ncol=5, nrow=length(unique(x.sp@cluster)) )
colnames(cluster_matrix) = c( "cluster_id", "WT1", "WT2", "KO1", "KO2" )
for ( i in 1:13 ) {
    cluster_matrix[i,1] = i
    cluster_matrix[i,2] = sum( x.sp@cluster == i & x.sp@sample == "WT1" )
    cluster_matrix[i,3] = sum( x.sp@cluster == i & x.sp@sample == "WT2" )
    cluster_matrix[i,4] = sum( x.sp@cluster == i & x.sp@sample == "KO1" )
    cluster_matrix[i,5] = sum( x.sp@cluster == i & x.sp@sample == "KO2" )
}

write.table( cluster_matrix, "analysis/2021-04-03_cluster_composition/cluster_matrix_2021-04-03.txt", sep="\t" )
