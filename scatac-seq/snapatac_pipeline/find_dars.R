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
load( "Rdata/call_ensemble_peaks_out.Rdata" )

#add peak matrix to snap object in R
x.sp = addPmatToSnap( x.sp )

#find differentially accessible peaks for all clusters

# pdf( paste0( "figures/",args[2],"_",args[3],"_DAR_peaks_viz_percluster.pdf"), width=8, height=5 )
# for ( i in levels( x.sp@cluster ) ) {
# DARs = findDAR(
#     obj=x.sp,
#     input.mat="pmat",
#     cluster.pos=i,
#     cluster.neg.method="knn",
#     test.method="exactTest",
#     bcv=0.1, #0.4 for human, 0.1 for mouse
#     seed.use=10
# );

# DARs$FDR = p.adjust(DARs$PValue, method="BH");
# idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
# write_delim( DARs, paste0( "peaks/diffpeaks/diffpeaks_cluster_",i,"_",args[3],".txt" ), delim="\t", col_names=T )

# par(mfrow = c(1, 2));
# plot(DARs$logCPM, DARs$logFC, 
#     pch=19, cex=0.1, col="grey", 
#     ylab="logFC", xlab="logCPM",
#     main=i
#   );
# points(DARs$logCPM[idy], 
#     DARs$logFC[idy], 
#     pch=19, 
#     cex=0.5, 
#     col="red"
# );

# abline(h = 0, lwd=1, lty=2);

# covs = Matrix::rowSums(x.sp@pmat);
# vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
# vals.zscore = (vals - mean(vals)) / sd(vals);

# plotFeatureSingle(
#     obj=x.sp,
#     feature.value=vals.zscore,
#     method="umap", 
#     main=i,
#     point.size=0.1, 
#     point.shape=19, 
#     down.sample=20000,
#     quantiles=c(0.01, 0.99)
# );
# }
# dev.off()

#find DARs for all clusters using snapatac default script
idy.ls = lapply(levels(x.sp@cluster), function(cluster_i){
	DARs = findDAR(
		obj=x.sp,
		input.mat="pmat",
		cluster.pos=cluster_i,
		cluster.neg=NULL,
		cluster.neg.method="knn",
		bcv=0.1,
		test.method="exactTest",
		seed.use=10
		)
	DARs$FDR = p.adjust(DARs$PValue, method="BH")
	idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0)
	write_delim( DARs, paste0( "peaks/DAR/DAR_cluster_",cluster_i,"_",args[3],".txt" ), delim="\t", col_names=T )
	if((x=length(idy)) < 2000L){
			PValues = DARs$PValue
			PValues[DARs$logFC < 0] = 1
			idy = order(PValues, decreasing=FALSE)[1:2000]
			rm(PValues) # free memory
	}
	idy
})

names(idy.ls) = levels(x.sp@cluster)

pdf( paste0( "figures/",args[2],"_",args[3],"_DAR_peaks_viz_allclusters.pdf"), width=12, height=10 )
par(mfrow = c(3, 4))
for(cluster_i in levels(x.sp@cluster)){
	print(cluster_i)
	idy = idy.ls[[as.integer(cluster_i)]]
	covs = Matrix::rowSums(x.sp@pmat);
	vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs
	vals.zscore = (vals - mean(vals)) / sd(vals)
	plotFeatureSingle(
		obj=x.sp,
		feature.value=vals.zscore,
		method="umap", 
		main=cluster_i,
		point.size=0.1, 
		point.shape=19, 
		down.sample=20000,
		quantiles=c(0.01, 0.99)
	)
}
dev.off()

print("snapatac pipeline: .. saving work ... ")

save.image( "Rdata/find_dars_out.Rdata" )