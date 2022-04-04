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
library("TFBSTools")
library("JASPAR2016")
library("JASPAR2018")
library("BiocParallel")
library("BuenColors")
library("cowplot")
library("pheatmap")
library("ggrepel")
library("scales")
library("RColorBrewer")
library("wesanderson")

load( "Rdata/run_chromvar_unwrapped_out.Rdata" )
# load( "Rdata/subclustered_postaire_out.Rdata" )

args = c( getwd(), "WT1_WT2_KO1_KO2_merged", "2021-06-12", "8" )

custom_colormaps = jdb_color_maps
names(custom_colormaps) = seq(1:17)

#use chromvar for motif analysis percell -- need to use snapatac_env_plusChromVar for addt'l packages
register(MulticoreParam(args[4]))

x.sp = makeBinary(x.sp, "pmat")

#####
#core chromvar method 
data.use = x.sp@pmat
peak.use = x.sp@peak
min.count = 10

idy = which( Matrix::colSums(data.use) >= min.count )
data.use = data.use[,idy,dropping=T]

peak.use = peak.use[idy]

rse <- SummarizedExperiment(
    assays = list( counts=t(data.use) ),
    rowRanges = peak.use,
    colData = DataFrame( Cell_Type=1:nrow(data.use), depth = Matrix::rowSums(data.use) )
)

#compute gc bias
rse <- addGCBias( rse, genome=BSgenome.Mmusculus.UCSC.mm10 )

#getting jaspar motifs individually
# motifs <- getJasparMotifs( collection="CORE", species="Mus musculus" )
# motifs <- getJasparMotifs( collection="PBM", species="Mus musculus" )
# motifs <- getJasparMotifs( collection="CORE", species="Homo sapiens" )
# motifs <- getJasparMotifs( collection="PBM", species="Mus musculus" )

#getting all jaspar motifs in one search
opts <- list()
opts[["species"]] <- 9606
opts[["collection"]] <- "CORE"
homoMatrix <- getMatrixSet( JASPAR2018, opts)
opts[["species"]] <- 10090
opts[["collection"]] <- "CORE"
musMatrix <- getMatrixSet( JASPAR2018, opts)
opts[["species"]] <- 10090
opts[["collection"]] <- "PBM"
musPBMMatrix <- getMatrixSet( JASPAR2018, opts)
motifs <- c( musMatrix, musPBMMatrix, homoMatrix )
if (!isTRUE(all.equal(TFBSTools::name(motifs), names(motifs)))) 
names(motifs) <- paste(names(motifs), TFBSTools::name(motifs), sep = "_")

#scanning motifs in the peaks
motif_mm <- matchMotifs( motifs, rse, genome=BSgenome.Mmusculus.UCSC.mm10 )

#calculate motif variability b/t cells
dev <- computeDeviations( object=rse, annotations=motif_mm )
dev_mat = t( assay(dev) )

#compute motif variability
dev <- dev[!duplicated( rownames(dev) ),]
variability <- computeVariability(dev)

##### visualization
#visualize results in boxplot and in umap
wt_idx <- x.sp@sample %in% c("WT1","WT2")
ko_idx <- x.sp@sample %in% c("KO1","KO2")

x.sp@metaData$cluster <- x.sp@cluster

for (motif_i in colnames(dev_mat) ) {
    # motif_i = "MA0080.4_SPI1"; #choose any given motif to viz
    dat = data.frame(x=factor(x.sp@metaData$cluster, levels=c( 6,11,12,1,7,8,9,10,4,5,13,3,2 ) ), y=dev_mat[,motif_i])
    y_quant = quantile( dat[,2], 1 )
    p <- ggplot2::ggplot(dat, aes(x=x, y=y, fill=x)) + 
        ggplot2::geom_hline( yintercept=0, linetype=3, color="gray" ) +
        ggplot2::geom_jitter( alpha=0.2, width=0.25, color="gray" ) +
        ggplot2::theme_classic() +
	    # ggplot2::geom_violin( draw_quantiles=0.5, width=0.5 ) + 
  	    ggplot2::geom_boxplot( width=0.35, outlier.shape=NA ) + 
	    ggplot2::xlab("cluster") +
	    ggplot2::ylab( "Relative enrichment\n(unbiased Z-score)" ) + 
        ggplot2::ylim( -y_quant,y_quant ) +
        ggplot2::scale_fill_manual( values = custom_colormaps ) +
        ggplot2::ggtitle( paste0( motif_i, " all" ) ) +
        ggplot2::theme( legend.position = "none" )
    q <- ggplot2::ggplot(dat[wt_idx,], aes(x=x, y=y, fill=x)) + 
        ggplot2::geom_hline( yintercept=0, linetype=3, color="gray" ) +
        ggplot2::geom_jitter( alpha=0.2, width=0.25, color="gray" ) +
        ggplot2::theme_classic() +
        # ggplot2::geom_violin( draw_quantiles=0.5, width=0.5 ) + 
        ggplot2::geom_boxplot( width=0.35, outlier.shape=NA ) + 
        ggplot2::xlab("cluster") +
        ggplot2::ylab( "Relative enrichment\n(unbiased Z-score)" ) + 
        ggplot2::ylim( -y_quant,y_quant ) +
        ggplot2::scale_fill_manual( values = custom_colormaps ) +
        ggplot2::ggtitle(paste0( motif_i, " WT" ) ) +
        ggplot2::theme( legend.position = "none" )
    r <- ggplot2::ggplot(dat[ko_idx,], aes(x=x, y=y, fill=x)) + 
        ggplot2::geom_hline( yintercept=0, linetype=3, color="gray" ) +
        ggplot2::geom_jitter( alpha=0.2, width=0.25, color="gray" ) +
        ggplot2::theme_classic() +
        # ggplot2::geom_violin( draw_quantiles=0.5, width=0.5 ) + 
        ggplot2::geom_boxplot( width=0.35, outlier.shape=NA ) + 
        ggplot2::xlab("cluster") +
        ggplot2::ylab( "Relative enrichment\n(unbiased Z-score)" ) + 
        ggplot2::ylim( -y_quant,y_quant ) +
        ggplot2::scale_fill_manual( values = custom_colormaps ) +
        ggplot2::ggtitle(paste0( motif_i, " KO" ) ) +
        ggplot2::theme(	legend.position = "none" )
    pdf( paste0( "figures/2020-10-01_chromvar/boxplot/",args[2],"_chromvar_",motif_i,"_",args[3],"_unwrapped.pdf"), height=7, width=8 )
    print( plot_grid(p,q,r, ncol=1) )
    dev.off()
    s <- ggplot2::ggplot( as.data.frame(x.sp@umap), ggplot2::aes( x=x.sp@umap[order(dev_mat[,motif_i]),1], y=x.sp@umap[order(dev_mat[,motif_i]),2] ) ) +
            ggplot2::geom_point( aes( color = dev_mat[order(dev_mat[,motif_i]), motif_i ],
                                      alpha = dev_mat[order(dev_mat[,motif_i]), motif_i ] ) ) +
            ggplot2::theme_classic() +
            ggplot2::xlab("UMAP-1") +
            ggplot2::ylab("UMAP-2") +
            ggplot2::scale_color_gradient( low="lightgoldenrodyellow", high="red",
                                            limits=c( 0.05*max(dev_mat[,motif_i]), 
                                                     0.95*max(dev_mat[,motif_i]) ),
                                            oob=squish ) +
            ggplot2::labs(col="gray") +
            ggplot2::ggtitle(motif_i) +
            ggplot2::theme(legend.position="none")
    pdf( paste0( "figures/2020-10-01_chromvar/umap/", args[2], "_", args[3],"_chromvar_unwrapped_",motif_i, "_umap_YlRd.pdf" ), height=3.2, width=3)
    print(s)
    dev.off()
}

#investigate most variable motifs
temp <- variability[,2] > 2
label_index <- ifelse( temp==T, as.vector( variability[,1] ), "" )
label_index <- ifelse( duplicated(label_index)==T, "", label_index)

temp <- variability[,1] %in% c("FOS","JUNB","JUND","RELA","Pou2f3","TP63","SNAI2","Hnf4a","HNF4G","GRHL1","GRHL2","FOXA1","Foxa2","Foxj3","SOX9","SPI1","SPIB","IRF4","IRF8","Myog","Myod1","ASCL1","Ascl2","CEBPA","CEBPB")
label_index <- ifelse( temp==T, as.vector( variability[,1] ), "" )
label_index <- ifelse( duplicated(label_index)==T, "", label_index)

# label_index <- ifelse( temp==T, substr( as.vector( variability[,1] ), 1, nchar( as.vector( variability[,1] ) ) - 2 ), "" )
label_index[ label_index %in% c("FOSL2","FOS::JUND","FOSL1::JUN",
                                "FOSL1","FOSL1::JUND","FOSL1::JUNB",
                                "POU2F1","FOSL2::JUN","BATF::JUN",
                                "NFKB2","POU5F1B") ] = ""

pdf( paste0( "figures/", args[2], "_", args[3],"_chromvar_unwrapped_variability_ggplot.pdf" ), height=3, width=5)
ggplot( as.data.frame(variability[,1]), aes( x=seq(1,length(variability[,2])), 
                                             y=variability[order(variability[,2], decreasing=T),2] ) ) +
    geom_point( size=2, shape=21, stroke=0.5, fill="gray", color="black" ) +
    # geom_errorbar( aes( ymin=variability[order(variability[,2], decreasing=T),3], 
    #                     ymax=variability[order(variability[,2], decreasing=T),4] ) ) +
    # geom_hline( yintercept=1.96, color="red", linetype="dashed" ) +
    geom_text_repel( label=label_index[order(variability[,2], decreasing=T)], xlim=c(50,700), force=1, segment.color="black", segment.alpha=0.25 ) +
    theme_classic() +
    xlab("JASPAR motifs") +
    ylab("Motif variability\n(sd of bias-corrected Z-score)")
dev.off()

#make heatmap of motifs

colnames(dev_mat)[ str_detect( colnames(dev_mat), "Foxa3" ) ]

features <- c(
  "MA0525.2_TP63","MA0101.1_REL",
  "MA0114.3_Hnf4a","MA0484.1_HNF4G",
  "MA0081.1_SPIB","PB0074.1_Sox8_1",
  "MA0148.3_FOXA1","MA0047.2_Foxa2",
  "PB0016.1_Foxj1_1","MA0600.2_RFX2",
  "MA0647.1_GRHL1","MA1105.1_GRHL2",
  "MA0627.1_Pou2f3"
)

keep <- features %in% colnames(dev_mat)
features <- features[keep]

sample_n <- 40
mmat_heatmap_mat <- dev_mat[x.sp@cluster==6,][ sample( nrow(dev_mat[x.sp@cluster==6,]), size=sample_n, replace=F ), ]
cluster_idx <- rep( 6, times=sample_n )

for (i in c(1,4,3,5,13,10,2)) {
    temp <- dev_mat[x.sp@cluster==i,][ sample( nrow(dev_mat[x.sp@cluster==i,]), size=sample_n, replace=F ), ]
    mmat_heatmap_mat <- rbind( mmat_heatmap_mat, temp )
    temp <- rep( i, times=sample_n )
    cluster_idx <- c( cluster_idx, temp )
}

keep <- colnames( mmat_heatmap_mat ) %in% features

mmat_heatmap_mat <- t( mmat_heatmap_mat[,keep] )
colnames(mmat_heatmap_mat) <- seq(1:ncol(mmat_heatmap_mat))

cluster_idx <- factor( cluster_idx, levels=c(6,1,4,3,5,13,10,2) )

breakslist = seq(0,0.15,by=0.01)
anno_colors <- list( cluster_idx = jdb_color_maps[1:13] )
names(anno_colors$cluster_idx) <- seq(1:13)
# row_labels <- substr( colnames(mmat_heatmap_mat), 8, nchar(colnames(mmat_heatmap_mat )) )

x <- mmat_heatmap_mat - rowMeans(mmat_heatmap_mat)

pdf( paste0( "figures/",args[2],"_",args[3],"_motif_heatmap_perCluster.pdf"), height=3, width=8)
pheatmap( x[features,],
          cluster_rows=F,
          cluster_cols=F,
          annotation_col=as.data.frame(cluster_idx),
          annotation_names_col=F,
          annotation_legend=F,
          annotation_colors=anno_colors,
          show_rownames=T,
          show_colnames=F,
        #   gaps_col=c(seq(100,800,by=100)),
        #   gaps_row=c(1,2,4,6,9,11),
        #   labels_row=row_labels
          breaks=breakslist,
          color=colorRampPalette((brewer.pal(name="Reds", n=9)))(length(breakslist))
)
dev.off()

#make heatmap of motifs in all cell of each cluster

colnames(dev_mat)[ str_detect( colnames(dev_mat), "Foxa3" ) ]

features <- c(
  "MA0525.2_TP63","MA0101.1_REL",
  "MA0114.3_Hnf4a","MA0484.1_HNF4G",
  "MA0081.1_SPIB","PB0074.1_Sox8_1",
  "MA0148.3_FOXA1","MA0047.2_Foxa2",
  "PB0016.1_Foxj1_1","MA0600.2_RFX2",
  "MA0647.1_GRHL1","MA1105.1_GRHL2",
  "MA0627.1_Pou2f3"
)

keep <- features %in% colnames(dev_mat)
features <- features[keep]

temp <- dev_mat[,features]
cluster_idx <- c(6,1,4,3,5,13,10,2)
mmat_heatmap_mat <- matrix( 0, nrow=length(cluster_idx), ncol=length(features) )
rownames(mmat_heatmap_mat) <- cluster_idx
colnames(mmat_heatmap_mat) <- features
count <- 0
for ( i in cluster_idx ) {
    count <- count+1
    print( paste0( "analyzing cluster ", i ) )
    print( paste0( "summing these motifs: ", features ) )
    mmat_heatmap_mat[count,] <- colMeans( temp[x.sp@cluster==i,] )
}

cluster_idx <- factor( cluster_idx, levels=c(6,1,4,3,5,13,10,2) )
mmat_heatmap_mat <- t( mmat_heatmap_mat )

breakslist = seq(0,0.15,by=0.01)
anno_colors <- list( cluster_idx = jdb_color_maps[1:13] )
names(anno_colors$cluster_idx) <- seq(1:13)
# row_labels <- substr( colnames(mmat_heatmap_mat), 8, nchar(colnames(mmat_heatmap_mat )) )

pdf( paste0( "figures/",args[2],"_",args[3],"_motif_heatmap_avgPerCluster_Reds.pdf"), height=3, width=8)
pheatmap( mmat_heatmap_mat,
          cluster_rows=F,
          cluster_cols=F,
          annotation_col=as.data.frame(cluster_idx),
          annotation_names_col=F,
          annotation_legend=F,
          annotation_colors=anno_colors,
          show_rownames=T,
          show_colnames=F,
          gaps_col=c(seq(1,length(cluster_idx)-1)),
        #   gaps_row=c(1,2,4,6,10,12),
        #   labels_row=row_labels,
          border_color=NA,
          breaks=breakslist,
          color=colorRampPalette((brewer.pal(name="Reds", n=9)))(length(breakslist))
        #   color=colorRampPalette(rev(brewer.pal(name="PRGn", n=9)))(length(breakslist))
        
)
dev.off()

#view individual motifs
features <- c(
  "MA0525.2_TP63","MA0101.1_REL",
  "MA0114.3_Hnf4a","MA0484.1_HNF4G",
  "MA0081.1_SPIB","PB0074.1_Sox8_1",
  "MA0148.3_FOXA1","MA0047.2_Foxa2",
  "PB0016.1_Foxj1_1","MA0600.2_RFX2",
  "MA0647.1_GRHL1","MA1105.1_GRHL2",
  "MA0627.1_Pou2f3"
)

wt_idx <- x.sp@sample %in% c("WT1","WT2") & x.sp@cluster %in% c(1)
ko_idx <- x.sp@sample %in% c("KO1","KO2") & x.sp@cluster %in% c(7)
keep <- wt_idx==T | ko_idx ==T
total_idx <- wt_idx + ko_idx*2
total_idx <- total_idx[keep]

for( anno_i in features ){
dat = data.frame(x=factor(x.sp@metaData$cluster, levels=c( 6,11,12,1,7,8,9,10,4,5,3,2 ) ), y=dev_mat[,anno_i]);
y_quant = quantile( dat[keep,2], 1 )
    s <- ggplot2::ggplot(dat[keep,], aes(x=x, y=y, fill=x)) + 
        ggplot2::geom_hline( yintercept=0, linetype=3, color="gray" ) +
        ggplot2::geom_jitter( alpha=0.2, width=0.25, color="gray" ) +
        ggplot2::theme_classic() +
        # ggplot2::geom_violin( draw_quantiles=0.5, width=0.5 ) + 
        ggplot2::geom_boxplot( width=0.35, outlier.shape=NA ) + 
        ggplot2::xlab("cluster") +
        ggplot2::ylab( "Relative enrichment\n(unbiased Z-score)" ) + 
        ggplot2::ylim( -y_quant,y_quant ) +
        ggplot2::scale_fill_manual( values = custom_colormaps ) +
        ggplot2::ggtitle(paste0( anno_i, "WT vs KO" ) ) +
        ggplot2::theme(	legend.position = "none" )
    pdf( paste0( "figures/",args[2],"_chromvar_",anno_i,"_",args[3],"_WT1vsKO7.pdf"), height=3, width=2 )
    print( s )
    dev.off()
}



#use annotations from bed
my_annotation_files <- c("./reference/aire_peaks.bed",
                         "reference/mm10.AIG.genes.bed",
                         "reference/mm10.ANG.genes.bed", 
                        #  "reference/WT_KO_CTCF_SMC1_MED1_Aire.peaks_with_adapter_IgG_annotated_2_d_given_CTCF_SMC1_common.bed", 
                        #  "reference/WT_H3K27ac_2_MSA_1_merge.peaks_with_adapter_IgG_super_enhancer_rose.bed",
                        #  "reference/WT_CTCF.peaks_with_adapter_IgG_factor.bed",
                         "/n/groups/cbdm-db/dam41/GSE114713_hollander/peaks/macs/2020-06-10/SRR7190777_h3k9me3_peaks.broadPeak",
                        #  "/n/groups/cbdm-db/dam41/GSE114713_hollander/peaks/macs/2020-06-10/SRR7190774_h3k9ac_peaks.narrowPeak",
                         "/n/groups/cbdm-db/dam41/GSE114713_hollander/peaks/macs/2020-06-10/SRR7190772_h3k4me1_peaks.narrowPeak",
                        #  "/n/groups/cbdm-db/dam41/GSE114713_hollander/peaks/macs/2020-06-10/SRR7190770_h3k4ac_peaks.narrowPeak",
                         "/n/groups/cbdm-db/dam41/GSE114713_hollander/peaks/macs/2020-06-10/SRR7190768_h3k27me3_peaks.broadPeak",
                         "/n/groups/cbdm-db/dam41/GSE53109_hollander/peaks/macs/2020-06-10/SRR1045005_h3k4me3_peaks.narrowPeak"
                        #  "/n/groups/cbdm-db/dam41/mec_scatac/reference/scratch/aire_h3k4me3_overlapping.bed",
                        #  "/n/groups/cbdm-db/dam41/mec_scatac/reference/scratch/aire_h3k4me1_overlapping.bed",
                        #  "/n/groups/cbdm-db/dam41/mec_scatac/reference/scratch/aire_h3k27me3_overlapping.bed" 
                        )

anno_ix <- getAnnotations(my_annotation_files, 
                          rowRanges = rowRanges(rse))

dev <- computeDeviations( object=rse, annotations=anno_ix )

rownames(dev) <- c("Aire","AIG","ANG","H3K9me3","H3K4me1","H3K27me3","H3K4me3")
dev_mat = t( assay(dev) )
colnames(dev_mat) <- c("Aire","AIG","ANG","H3K9me3","H3K4me1","H3K27me3","H3K4me3")

wt_idx <- x.sp@sample %in% c("WT1","WT2") & x.sp@cluster %in% c(1)
ko_idx <- x.sp@sample %in% c("KO1","KO2") & x.sp@cluster %in% c(7)
keep <- wt_idx==T | ko_idx ==T
total_idx <- wt_idx + ko_idx*2
total_idx <- total_idx[keep]

print("Aire")
test <- wilcox.test( dev_mat[wt_idx,"Aire"],dev_mat[ko_idx,"Aire"],alternative="two.sided",paired=F )$p.value / 7
print(test)
print("AIG")
test <- wilcox.test( dev_mat[wt_idx,"AIG"],dev_mat[ko_idx,"AIG"],alternative="two.sided",paired=F )$p.value / 7
print(test)
print("ANG")
test <- wilcox.test( dev_mat[wt_idx,"ANG"],dev_mat[ko_idx,"ANG"],alternative="two.sided",paired=F )$p.value / 7
print(test)
print("H3K9me3")
test <- wilcox.test( dev_mat[wt_idx,"H3K9me3"],dev_mat[ko_idx,"H3K9me3"],alternative="two.sided",paired=F )$p.value / 7
print(test)
print("H3K4me1")
test <- wilcox.test( dev_mat[wt_idx,"H3K4me1"],dev_mat[ko_idx,"H3K4me1"],alternative="two.sided",paired=F )$p.value / 7
print(test)
print("H3K27me3")
test <- wilcox.test( dev_mat[wt_idx,"H3K27me3"],dev_mat[ko_idx,"H3K27me3"],alternative="two.sided",paired=F )$p.value / 7
print(test)
print("H3K4me3")
test <- wilcox.test( dev_mat[wt_idx,"H3K4me3"],dev_mat[ko_idx,"H3K4me3"],alternative="two.sided",paired=F )$p.value / 7
print(test)

for( anno_i in colnames(dev_mat) ){
dat = data.frame(x=factor(x.sp@metaData$cluster, levels=c( 6,11,12,1,7,8,9,10,4,5,3,2 ) ), y=dev_mat[,anno_i]);
y_quant = quantile( dat[,2], 0.99 )
    s <- ggplot2::ggplot(dat[keep,], aes(x=x, y=y, fill=x)) + 
        ggplot2::geom_hline( yintercept=0, linetype=3, color="gray" ) +
        ggplot2::geom_jitter( alpha=0.2, width=0.25, color="gray" ) +
        ggplot2::theme_classic() +
        # ggplot2::geom_violin( draw_quantiles=0.5, width=0.5 ) + 
        ggplot2::geom_boxplot( width=0.35, outlier.shape=NA ) + 
        ggplot2::xlab("cluster") +
        ggplot2::ylab( "Relative enrichment\n(unbiased Z-score)" ) + 
        ggplot2::ylim( -y_quant,y_quant ) +
        ggplot2::scale_fill_manual( values = custom_colormaps ) +
        ggplot2::ggtitle(paste0( anno_i, "WT vs KO" ) ) +
        ggplot2::theme(	legend.position = "none" )
    pdf( paste0( "figures/",args[2],"_chromvar_",anno_i,"_",args[3],"_WT1vsKO7.pdf"), height=3, width=2 )
    print( s )
    dev.off()
}

##save work
# save.image( "Rdata/run_chromvar_unwrapped_out.Rdata" )