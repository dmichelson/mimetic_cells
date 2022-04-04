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
suppressMessages( suppressWarnings( library("JASPAR2018") ) )
suppressMessages( suppressWarnings( library("BiocParallel") ) )

print( "snapatac pipeline: .. loading data ... ")
load( "Rdata/find_dars_out.Rdata" )

args = commandArgs(trailingOnly=TRUE)

#use chromvar for motif analysis percell -- need to use snapatac_env_plusChromVar for addt'l packages
register(MulticoreParam(args[4]))

print( "snapatac pipeline: .. making binary pmat ... " )
x.sp = makeBinary(x.sp, "pmat");

print( "snapatac pipeline: .. running chromVAR ... " )
x.sp@mmat = runChromVAR(
    obj=x.sp,
    input.mat="pmat",
    genome=BSgenome.Mmusculus.UCSC.mm10,
    min.count=10,
    species="Mus musculus"
)
x.sp;

print( "snapatac pipeline: .. visualizing motifs ... " )
for (motif_i in colnames(x.sp@mmat)) {
# motif_i = "MA1684.1_Foxn1"; #choose any given motif to viz
dat = data.frame(x=x.sp@metaData$cluster, y=x.sp@mmat[,motif_i]);
p <- ggplot2::ggplot(dat, aes(x=x, y=y, fill=x)) + 
    ggplot2::theme_classic() +
	ggplot2::geom_violin() + 
	ggplot2::xlab("cluster") +
	ggplot2::ylab("motif enrichment") + 
	ggplot2::ggtitle(motif_i) +
	ggplot2::theme(
		  plot.margin = margin(5,1,5,1, "cm"),
		  axis.text.x = element_text(angle = 90, hjust = 1),
		  axis.ticks.x=element_blank(),
		  legend.position = "none"
)

pdf( paste0( "figures/chromvar/",args[2],"_chromvar_",motif_i,"_",args[3],".pdf"), height=6, width=10 )
print(p)
dev.off()
}

print( "snapatac pipeline: .. saving work ... " )

save.image( "Rdata/run_chromvar_out.Rdata" )