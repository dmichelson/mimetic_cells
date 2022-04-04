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
suppressMessages( suppressWarnings( library("Seurat") ) )
suppressMessages( suppressWarnings( library("chromVAR") ) )
suppressMessages( suppressWarnings( library("motifmatchr") ) )
suppressMessages( suppressWarnings( library("SummarizedExperiment") ) )
suppressMessages( suppressWarnings( library("BSgenome.Mmusculus.UCSC.mm10") ) )
suppressMessages( suppressWarnings( library("JASPAR2016") ) )
suppressMessages( suppressWarnings( library("BiocParallel") ) )

print( "snapatac pipeline: .. loading data ... ")
load( "Rdata/cluster_viz_subset_out.Rdata" )

#integrate scRNA-seq
mec.rna = readRDS("scrna/mec_seurat_table_2020-04-03.rds")
mec.rna$tech = "rna"
variable.genes = VariableFeatures(object = mec.rna)
genes.df = read.table("reference/gencode.vM16.gene.bed")
genes.gr = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4])
genes.sel.gr = genes.gr[which(genes.gr$name %in% variable.genes)]

## reload the bmat, this is optional but recommended
print( "snapatac pipeline: .. reloading bmat ..." )
x.sp = addBmatToSnap(x.sp);
x.sp = createGmatFromMat(
    obj=x.sp, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores= as.integer( args[4] )
)

#proceed with scrna integration
print( "snapatac pipeline: .. integrating scrna and scatac ..." )
mec.atac <- snapToSeurat(
    obj=x.sp, 
    eigs.dims=1:20, 
    norm=TRUE,
    scale=TRUE
)
transfer.anchors <- FindTransferAnchors(
    reference = mec.rna, 
    query = mec.atac, 
    features = variable.genes, 
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    reduction = "cca"
)
celltype.predictions <- TransferData(
    anchorset = transfer.anchors, 
    refdata = mec.rna@active.ident,
    weight.reduction = mec.atac[["SnapATAC"]],
    dims = 1:20
)
x.sp@metaData$predicted.id = celltype.predictions$predicted.id;
x.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max);

#create pseudo multiomics cells

print( "snapatac pipeline: .. adding pseudo RNA data to scatac dataset ..." )

refdata <- GetAssayData(
    object = mec.rna, 
    assay = "RNA", 
    slot = "data"
)
imputation <- TransferData(
    anchorset = transfer.anchors, 
    refdata = refdata, 
    weight.reduction = mec.atac[["SnapATAC"]], 
    dims = 1:20
)
x.sp@gmat = t(imputation@data)
rm(imputation); # free memory
rm(refdata);    # free memory
rm(mec.rna);   # free memory
rm(mec.atac); # free memory

#remove cells with a low prediction score
print( "snapatac pipeline: .. plotting prediction score.. identifying low prediction score cells ..." )
pdf( paste0("figures/", args[2], "_scrna_prediction_score_", args[3], ".pdf" ), width=5, height=5 )
hist(
    x.sp@metaData$predict.max.score, 
    xlab="prediction score", 
    col="lightblue", 
    xlim=c(0, 1),
    main="MEC scATAC integration w scRNA, prediction score"
)
abline(v=0.2, col="red", lwd=2, lty=2)
table(x.sp@metaData$predict.max.score > 0.2)
dev.off()

x.sp[x.sp@metaData$predict.max.score > 0.2,]
# x.sp = x.sp[x.sp@metaData$predict.max.score > 0.2,]

pdf( paste0("figures/", args[2],"_", args[3], "_scrna_scatac_predicted_id.pdf" ), width=5, height=5 )
plotViz(
    obj=x.sp,
    method="umap", 
    main="MEC scATAC integration w scRNA, UMAP",
    point.color=x.sp@metaData[,"predicted.id"], 
    point.size=0.5, 
    point.shape=19, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    down.sample=20000,
    legend.add=FALSE
)
dev.off()

#viz marker genes from big gmat

marker.genes = c(
    "H2-Aa", "Cd74", "Aire", "Cd80",
    "Pigr", "Ly6a","Pdpn","Ctla4","Il10"
);


pdf(paste( "figures/", args[2], "_", args[3], "_cluster_viz_marker_genes.pdf", sep="" ), height=13, width=10)
par(mfrow = c(3, 3))
for(i in 1:9){
    plotFeatureSingle(
        obj=x.sp,
        feature.value=x.sp@gmat[, marker.genes[i]],
        method="umap", 
        main=marker.genes[i],
        point.size=0.1, 
        point.shape=19, 
        down.sample=10000,
        quantiles=c(0, 1)
	)
}
dev.off()

print( "snapatac pipeline: .. saving work ..." )

save.image( "Rdata/integrate_scrna_out.Rdata" )
