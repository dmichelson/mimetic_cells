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

#set sample ids, files, barcodes
sample_list = c(
    "WT_Aire_MEChi1_10xATAC_5K_REA09522",
    "WT2_reprep_Aire_MEChi1_10xATAC_5K",
    "KO_Aire_MEChi1_10xATAC_5K_REA09522",
    "KO2_reprep_Aire_MEChi1_10xATAC_5K"
)

snap.files = c(
    "/n/groups/cbdm-db/dam41/mec_scatac/alignments/WT_Aire_MEChi1_10xATAC_5K_REA09522.snap",
	"/n/groups/cbdm-db/dam41/mec_scatac/alignments/WT2_reprep_Aire_MEChi1_10xATAC_5K.snap",
	"/n/groups/cbdm-db/dam41/mec_scatac/alignments/KO_Aire_MEChi1_10xATAC_5K_REA09522.snap",
	"/n/groups/cbdm-db/dam41/mec_scatac/alignments/KO2_reprep_Aire_MEChi1_10xATAC_5K.snap"
)
	
sample.names = c(
    "WT1",
    "WT2",
	"KO1",
	"KO2"
)
	
barcode.files = c(
    "/n/groups/cbdm-db/dam41/mec_scatac/data/200305_A00794_0197_BHHM5WDRXX/count/WT_Aire_MEChi1_10xATAC_5K_REA09522/outs/singlecell.csv",
    "/n/groups/cbdm-db/dam41/mec_scatac/data/200305_A00794_0197_BHHM5WDRXX/count/WT2_reprep_Aire_MEChi1_10xATAC_5K/outs/singlecell.csv",
    "/n/groups/cbdm-db/dam41/mec_scatac/data/200305_A00794_0197_BHHM5WDRXX/count/KO_Aire_MEChi1_10xATAC_5K_REA09522/outs/singlecell.csv",
    "/n/groups/cbdm-db/dam41/mec_scatac/data/200305_A00794_0197_BHHM5WDRXX/count/KO2_reprep_Aire_MEChi1_10xATAC_5K/outs/singlecell.csv"
)

#create snap file	
x.sp.ls = lapply(seq(snap.files), function(i){
    createSnap(
        file=snap.files[i],
        sample=sample.names[i]
    );
})

names(x.sp.ls) = sample.names;

#read in barcodes, cleanup, calculate stats
barcode.ls = lapply(seq(snap.files), function(i){
    barcodes = read.csv(
        barcode.files[i], 
        head=TRUE
    );
    barcodes = barcodes[2:nrow(barcodes),];
	barcodes$barcode = substr(as.character(barcodes$barcode),1,nchar(as.character(barcodes$barcode))-2) #trim barcode in metadata to match my snap file
    barcodes$logUMI = log10(barcodes$passed_filters + 1);
    barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
    barcodes
})

#view barcode stats (frip, umi)
plots = lapply(seq(snap.files), function(i){
    p1 = ggplot(
        barcode.ls[[i]], 
        aes(x=logUMI, y=promoter_ratio)) + 
        geom_point(size=0.3, col="grey") +
        theme_classic()	+
        ggtitle(sample.names[[i]]) +
        ylim(0, 1) + xlim(0, 6) + 
        labs(x = "log10(UMI)", y="promoter ratio") +
        geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
        geom_hline(yintercept=0.75, linetype="dashed", color = "red") +
        geom_vline(xintercept=3.5, linetype="dashed", color = "red") +
        geom_vline(xintercept=5, linetype="dashed", color = "red")
        p1
})

pdf( paste("figures/",args[2],"_",args[3],"_frip_vs_umi.pdf", sep=""), width=8, height=8 )
plots
dev.off()

#select which barcodes to keep
cutoff.logUMI.low = c(3.5,3.5,3.5,3.5);
cutoff.logUMI.high = c(5,5,5,5);
cutoff.FRIP.low = c(0.2,0.2,0.2,0.2);
cutoff.FRIP.high = c(0.75,0.75,0.75,0.75);

barcode.ls = lapply(seq(snap.files), function(i){
    barcodes = barcode.ls[[i]];
    idx = which(
        barcodes$logUMI >= cutoff.logUMI.low[i] & 
        barcodes$logUMI <= cutoff.logUMI.high[i] & 
        barcodes$promoter_ratio >= cutoff.FRIP.low[i] &
        barcodes$promoter_ratio <= cutoff.FRIP.high[i]
    );
    barcodes[idx,]
});

x.sp.ls = lapply(seq(snap.files), function(i){
    barcodes = barcode.ls[[i]];
    x.sp = x.sp.ls[[i]];
    barcode.shared = intersect(x.sp@barcode, barcodes$barcode);
    x.sp = x.sp[match(barcode.shared, x.sp@barcode),];
    barcodes = barcodes[match(barcode.shared, barcodes$barcode),];
    x.sp@metaData = barcodes;
    x.sp
})

names(x.sp.ls) = sample.names;

x.sp.ls

#reduce snaps to single snap file with sample ID as metadata
x.sp = Reduce(snapRbind, x.sp.ls);
x.sp@metaData["sample"] = x.sp@sample;
x.sp

#write barcodes to file
count = 1
for (i in sample.names ) {
    j = count
    idx = x.sp@sample==i
    write_barcodes = as.data.frame(x.sp@barcode[idx])
    colnames(write_barcodes) = "barcode"
    write_delim(write_barcodes, paste(args[1],"/alignments/barcodes/",i,"_barcodes_umi",cutoff.logUMI.low[j],"-",cutoff.logUMI.high[j],"_frip",cutoff.FRIP.low[j],"-",cutoff.FRIP.high[j],"_",args[3],".txt",sep=""), 
    delim="\t", col_names=T )
    rm(idx)
    rm(write_barcodes)
    count = count + 1
}

#give barcode stats
sink( paste("alignments/qc/barcodes_passing_filtering_",args[3],".txt",sep="") )
for (i in sample.names ) {
    print( i )
    print( x.sp[x.sp@sample==i] )
    print( paste0("median logUMI=", median( x.sp@metaData[x.sp@sample == i, "logUMI" ] ) ) )
    print( paste0("median promoter_ratio=", median( x.sp@metaData[x.sp@sample == i, "promoter_ratio" ] ) ) )
    print( paste0("mean logUMI=", mean( x.sp@metaData[x.sp@sample == i, "logUMI" ] ) ) )
    print( paste0("mean promoter_ratio=", mean( x.sp@metaData[x.sp@sample == i, "promoter_ratio" ] ) ) )
    print( paste0("sd logUMI=", sd( x.sp@metaData[x.sp@sample == i, "logUMI" ] ) ) )
    print( paste0("sd promoter_ratio=", sd( x.sp@metaData[x.sp@sample == i, "promoter_ratio" ] ) ) )
}
sink()

print( "snapatac pipeline: .. saving work ... " )

save.image( "Rdata/preprocess_snap_out.Rdata" )