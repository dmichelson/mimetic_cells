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

temp <- read.delim("alignments/barcodes/WT1_barcodes_post_cluster_removal_2020-04-03.txt", sep="\t", header=T)
temp_table <- cbind( rep( "WT1",times=nrow(temp) ), temp )
colnames(temp_table) <- c("sample", "barcode")
temp_full_table <- temp_table
temp <- read.delim("alignments/barcodes/WT2_barcodes_post_cluster_removal_2020-04-03.txt", sep="\t", header=T)
temp_table <- cbind( rep( "WT2",times=nrow(temp) ), temp )
colnames(temp_table) <- c("sample", "barcode")
temp_full_table <- rbind(temp_full_table, temp_table)
temp <- read.delim("alignments/barcodes/KO1_barcodes_post_cluster_removal_2020-04-03.txt", sep="\t", header=T)
temp_table <- cbind( rep( "KO1",times=nrow(temp) ), temp )
colnames(temp_table) <- c("sample", "barcode")
temp_full_table <- rbind(temp_full_table, temp_table)
temp <- read.delim("alignments/barcodes/KO2_barcodes_post_cluster_removal_2020-04-03.txt", sep="\t", header=T)
temp_table <- cbind( rep( "KO2",times=nrow(temp) ), temp )
colnames(temp_table) <- c("sample", "barcode")
temp_full_table <- rbind(temp_full_table, temp_table)

idx <- logical( length=length(x.sp@barcode) )
for (i in sample.names) {
keep <- x.sp@barcode %in% temp_full_table[temp_full_table[,1]==i, 2]
idx <- idx + keep
}
idx <- idx > 0

#per barcode
pdf(paste( "figures/", args[2], "_", args[3], "_clusters_to_remove_viz.pdf", sep="" ), height=5, width=6)
ggplot2::ggplot( as.data.frame(x.sp@umap), ggplot2::aes( x=x.sp@umap[,1], y=x.sp@umap[,2] ) ) +
    ggplot2::geom_point( aes( color = idx , alpha=0.5 ) )+
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP-1") +
    ggplot2::ylab("UMAP-2") +
    ggplot2::scale_colour_manual( values=c("TRUE" = "blue", "FALSE" = "red") ) +
    ggplot2::theme( legend.position="none" ) +
    ggplot2::ggtitle("RED=REMOVE,BLUE=KEEP") +
    ggplot2::xlim(-15,11) +
    ggplot2::ylim(-14,7)
dev.off()

#remove per umap
idx = x.sp@umap[,1] < (-10) | ( x.sp@umap[,1] < (0) & x.sp@umap[,2] < (-7.5) ) | ( x.sp@umap[,1] > (7.5) & x.sp@umap[,2] < (-5) )
x.sp = x.sp[!idx,]

#remove per cluster
idx = x.sp@cluster %in% c(15,13,12)
x.sp = x.sp[!idx,]

rm(temp)
rm(temp_table)
rm(temp_full_table)

#viz cells post removal

pdf(paste( "figures/", args[2], "_", args[3], "_clusters_post_removal_viz.pdf", sep="" ), height=5, width=6)
ggplot2::ggplot( as.data.frame(x.sp@umap), ggplot2::aes( x=x.sp@umap[,1], y=x.sp@umap[,2] ) ) +
    ggplot2::geom_point( aes( color = "blue" , alpha=0.5 ) )+
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP-1") +
    ggplot2::ylab("UMAP-2") +
    ggplot2::theme( legend.position="none" ) +
    ggplot2::ggtitle("Clusters post-removal") +
    ggplot2::xlim(-15,11) +
    ggplot2::ylim(-14,7)
dev.off()

#write post cluster filter barcodes to file
count = 1
for (i in sample.names ) {
    j = count
    idx = x.sp@sample==i
    write_barcodes = as.data.frame(x.sp@barcode[idx])
    colnames(write_barcodes) = "barcode"
    write_delim(write_barcodes, paste0(args[1],"/alignments/barcodes/",i,"_barcodes_post_cluster_removal_",args[3],".txt"), 
    delim="\t", col_names=T )
    rm(idx)
    rm(write_barcodes)
    count = count + 1
}

save.image( "Rdata/remove_clusters_out.Rdata" )