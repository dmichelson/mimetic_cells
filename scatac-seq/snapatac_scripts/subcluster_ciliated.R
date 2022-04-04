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

load( "Rdata/call_ensemble_peaks_out.Rdata" )

args = c( getwd(), "WT1_WT2_KO1_KO2_merged", "2020-10-01", "8" )

custom_colormaps = jdb_color_maps
names(custom_colormaps) = seq(1:17)

#use snaptools built in dim reduction and clustering to examine similarity
feature_name = "bins"
samples = "WT_KO"
clusters = "5"
sample_id = c("WT1","WT2","KO1","KO2")
cluster = c("5")

idx = x.sp@sample %in% sample_id
sub.x.sp = x.sp[idx,,]
idx = sub.x.sp@cluster %in% cluster
sub.x.sp = sub.x.sp[idx,,]

sub.x.sp = makeBinary(sub.x.sp, mat="bmat")

sub.x.sp = runDimReduct(
    obj=sub.x.sp,
    input.mat="bmat", 
    pc.num=20,
    method="svd",
    center=T,
    scale=T,
    seed.use=10
)

pdf(paste0( "figures/", feature_name, "_", samples,"_", clusters, "_", args[3], "_dim_reduction_dimwise.pdf" ), height=4, width=5)
plotDimReductPW(
    obj=sub.x.sp, 
    eigs.dims=1:20,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
    )
dev.off()

sub.x.sp = runViz(
    obj=sub.x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:6, 
    method="umap",
    seed.use=10
)

k=2
kmeans_vector <- kmeans( sub.x.sp@umap, 
                         centers=k, 
                         iter.max=1000, 
                         algorithm="Hartigan-Wong" )

sub.x.sp@cluster <- as.factor(kmeans_vector$cluster)
sub.x.sp@metaData$cluster = sub.x.sp@cluster

pdf(paste0( "figures/", feature_name, "_", samples,"_", clusters, "_subsetted_", args[3],"_k-",k, ".pdf", sep="" ), height=5, width=5)
plotViz(
    obj=sub.x.sp,
    method="umap", 
    main=paste("MEC Cluster, k=",k,sep=""),
    point.color=sub.x.sp@cluster, 
    point.size=1, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=F,
    text.size=1.5,
    text.color="black",
    text.halo.add=F,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=20000,
    legend.add=FALSE
	)
dev.off()

gene="Foxj1"
pdf(paste0( "figures/", feature_name, "_", samples,"_", clusters, "_subsetted_", args[3],"_",gene, ".pdf", sep="" ), height=5, width=5)
plotFeatureSingle(
    obj=sub.x.sp,
    method="umap", 
    main=paste("MEC Cluster, ",gene,sep=""),
    feature.value=sub.x.sp@gmat[,gene], 
    point.size=0.5, 
    point.shape=19, 
    down.sample=20000,
	)
dev.off()


#integrate new clusters into full clustering
levels(x.sp@cluster) <- seq( 1, 13 )
idx <- x.sp@cluster==5 & x.sp@barcode %in% sub.x.sp@barcode
for (i in 1:sum(idx)) {
    if ( sub.x.sp@cluster[i]==1 ) {
        temp <- x.sp@barcode[idx]==sub.x.sp@barcode[i]
        x.sp@cluster[idx][temp] <- 13
    }
}

# combine clusters (optional)
# idx <- x.sp@cluster %in% c(7,8,9)
# for (i in 1:length(idx)) {
#     if ( idx[i]==T ) {
#         x.sp@cluster[i] <- 7
#     }
# }

# count = 0
# temp <- x.sp@cluster
# for (i in unique(x.sp@cluster)) {
#     count <- count+1
#     idx <- x.sp@cluster==i
#     temp[idx] <- count
# }
# x.sp@cluster <- temp

s <- ggplot2::ggplot( as.data.frame(x.sp@umap), ggplot2::aes( x=x.sp@umap[,1], y=x.sp@umap[,2] ) ) +
        ggplot2::geom_point( size=0.2, alpha=0.75, aes( color=x.sp@cluster ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        ggplot2::ggtitle("Cluster 5") +
        ggplot2::scale_color_manual( values=custom_colormaps ) +
        ggplot2::theme(legend.position="none")
pdf( paste0( "figures/", args[2], "_", args[3],"_umap_plus13.pdf" ), height=3.2, width=3)
print(s)
dev.off()

#run macs2
clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 0)]
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i])
    runMACS(
        obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
        output.prefix=paste0("peaks/",args[3],"/",args[2],"_",args[3],"_peaks", gsub(" ", "_", clusters.sel)[i]),
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

#make combined peak file
setwd( paste0("peaks/",args[3]) )
peaks.names = system("ls | grep narrowPeak", intern=TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
    peak.df = read.table(x)
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls))
peak.gr

setwd("..")
setwd("..")
#save combined peaks file
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = paste0("peaks/",args[3],"/",args[2],"_",args[3],"_peaks_combined.bed"),append=FALSE,
	quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
	row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
	fileEncoding = ""
)

#load pmat
x.sp = addPmatToSnap( x.sp )

feature_list <- peaks.df
feature_list.gr = GRanges( feature_list$seqnames, IRanges( feature_list$start, feature_list$end ) )
idy = queryHits(
    findOverlaps(x.sp@peak, feature_list.gr)
);

x.sp <- x.sp[,idy,"pmat"]

#find differentially accessible peaks for all clusters
pdf( paste0( "figures/",args[2],"_",args[3],"_DAR_peaks_viz_percluster.pdf"), width=8, height=5 )
for ( i in levels( x.sp@cluster ) ) {
DARs = findDAR(
    obj=x.sp,
    input.mat="pmat",
    cluster.pos=i,
    cluster.neg.method="knn",
    test.method="exactTest",
    bcv=0.1, #0.4 for human, 0.1 for mouse
    seed.use=10
);

DARs$FDR = p.adjust(DARs$PValue, method="BH");
idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0);
write_delim( DARs, paste0( "peaks/", args[3],"/diffpeaks_cluster_",i,"_",args[3],".txt" ), delim="\t", col_names=T )

par(mfrow = c(1, 2));
plot(DARs$logCPM, DARs$logFC, 
    pch=19, cex=0.1, col="grey", 
    ylab="logFC", xlab="logCPM",
    main=i
  );
points(DARs$logCPM[idy], 
    DARs$logFC[idy], 
    pch=19, 
    cex=0.5, 
    col="red"
);

abline(h = 0, lwd=1, lty=2);

covs = Matrix::rowSums(x.sp@pmat);
vals = Matrix::rowSums(x.sp@pmat[,idy]) / covs;
vals.zscore = (vals - mean(vals)) / sd(vals);

plotFeatureSingle(
    obj=x.sp,
    feature.value=vals.zscore,
    method="umap", 
    main=i,
    point.size=0.1, 
    point.shape=19, 
    down.sample=20000,
    quantiles=c(0.01, 0.99)
);
}
dev.off()

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
	if((x=length(idy)) < 2000L){
			PValues = DARs$PValue
			PValues[DARs$logFC < 0] = 1
			idy = order(PValues, decreasing=FALSE)[1:2000]
			rm(PValues) # free memory
	}
	idy
})

names(idy.ls) = levels(x.sp@cluster)

#find motifs in each cluster's DARs using homer
system("which findMotifsGenome.pl")
for ( i in levels(x.sp@cluster) ) {
    motifs = runHomer(
	    x.sp[,idy.ls[[as.integer(i)]],"pmat"], 
    	mat = "pmat",
    	path.to.homer = "/n/scratch3/users/d/dam41/homer/bin/findMotifsGenome.pl",
    	result.dir = paste0("peaks/",args[3],"/homer/cluster_",i),
    	num.cores=args[4],
    	genome = 'mm10',
    	motif.length = 10,
    	scan.size = 300,
    	optimize.count = 2,
    	background = 'automatic',
    	local.background = FALSE,
    	only.known = TRUE,
    	only.denovo = FALSE,
    	fdr.num = 5,
    	cache = 100,
    	overwrite = TRUE,
    	keep.minimal = FALSE,
    )
}

save.image ( "Rdata/subclustered_postaire_out.Rdata" )