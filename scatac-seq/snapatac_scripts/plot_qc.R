#!/usr/bin/env Rscript

#snapatac Rscript
#feed inputs as "Rscript run_snapatac_generic.R *wd *snap_file_name *date *n_cores [ *arg6...]"

#outside the environment, load homer and macs2 into path. Homer is in scratch directory and may need to be reinstalled...

#set up environment
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

load( "Rdata/cluster_viz_subset_out.Rdata" )

args =c( getwd(), "WT1_WT2_KO1_KO2_merged", "2020-05-31", "8" )

pdf("figures/read_depth_2020-05-31.pdf", height=5, width=6)	
ggplot2::ggplot( x.sp@metaData, aes( x=x.sp@umap[,1], y=x.sp@umap[,2] )  ) +
    ggplot2::geom_point( aes( color=log(x.sp@metaData[,"passed_filters"]+1,10) ) ) +
    viridis::scale_color_viridis( limits=c(3,5) ) +
    ggplot2::theme( legend.position="right" ) +
    ggplot2::labs( color="log(UMI)" ) +
    ggplot2::ggtitle( "Read depth" ) +
    ggplot2::theme_classic()
dev.off()

pdf("figures/duplicates_2020-05-31.pdf", height=5, width=6)	
ggplot2::ggplot( x.sp@metaData, aes( x=x.sp@umap[,1], y=x.sp@umap[,2] )  ) +
    ggplot2::geom_point( aes( color=x.sp@metaData$duplicate / x.sp@metaData$total ) ) +
    viridis::scale_color_viridis( limits=c(0,1) ) +
    ggplot2::theme( legend.position="right" ) +
    ggplot2::labs( color="fraction\nduplicates" ) +
    ggplot2::ggtitle( "Duplicates" ) +
    ggplot2::theme_classic()
dev.off()

pdf("figures/fripeak_2020-05-31.pdf", height=5, width=6)	
ggplot2::ggplot( x.sp@metaData, aes( x=x.sp@umap[,1], y=x.sp@umap[,2] )  ) +
    ggplot2::geom_point( aes( color=x.sp@metaData$peak_region_fragments / x.sp@metaData$passed_filters ) ) +
    viridis::scale_color_viridis( limits=c(0,1) ) +
    ggplot2::theme( legend.position="right" ) +
    ggplot2::labs( color="fraction\nin peaks" ) +
    ggplot2::ggtitle( "FRiP" ) +
    ggplot2::theme_classic()
dev.off()

pdf("figures/read_depth_2021-05-18.pdf", height=5, width=6)	
ggplot2::ggplot( x.sp@metaData  ) +
    ggplot2::geom_histogram( aes( x=log(x.sp@metaData[,"passed_filters"]+1,10) ), binwidth=0.2, color="black", fill="gray" ) +
    ggplot2::xlim(3,5) +
    ggplot2::theme( legend.position="none" ) +
    ggplot2::xlab( "log(UMI)" ) +
    ggplot2::ggtitle( "Read depth" ) +
    ggplot2::theme_bw()
dev.off()

pdf("figures/duplicates_2021-05-18.pdf", height=5, width=6)	
ggplot2::ggplot( x.sp@metaData  ) +
    ggplot2::geom_histogram( aes( x=x.sp@metaData$duplicate / x.sp@metaData$total ), binwidth=0.1, color="black", fill="gray" ) +
    ggplot2::xlim(0,1) +
    ggplot2::theme( legend.position="none" ) +
    ggplot2::xlab( "fraction\nduplicates" ) +
    ggplot2::ggtitle( "Duplicates" ) +
    ggplot2::theme_bw()
dev.off()

pdf("figures/fripeak_2021-05-18.pdf", height=5, width=6)	
ggplot2::ggplot( x.sp@metaData ) +
    ggplot2::geom_histogram( aes( x=x.sp@metaData$peak_region_fragments / x.sp@metaData$passed_filters ), binwidth=0.1, color="black", fill="gray" ) +
    ggplot2::xlim(0,1) +
    ggplot2::theme( legend.position="none" ) +
    ggplot2::xlab( "fraction\nin peaks" ) +
    ggplot2::ggtitle( "FRiP" ) +
    ggplot2::theme_bw()
dev.off()