library(tidyverse)
library(Seurat)
library(Matrix)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(future)
library(tidyverse)

# plan("multiprocess", workers=4)

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_100821")

date="20211124"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value

#####
#read gene expression data
dat <- Read10X( data.dir="data_cr/filtered_feature_bc_matrix" )
seurat_table <- CreateSeuratObject(dat, min.cells = 3, project="MEClo2", min.features=0)

rm(dat)

#read in hashes

hashtag_table <- readMM("data_cr/hash/umi_count/matrix.mtx.gz")
rownames(hashtag_table) <- read.table("data_cr/hash/umi_count/features.tsv.gz")$V1
colnames(hashtag_table) <- read.table("data_cr/hash/umi_count/barcodes.tsv.gz")$V1
colnames(hashtag_table) <- paste0(colnames(hashtag_table),"-1")

#integrate hashes
length(intersect(colnames(seurat_table), colnames(hashtag_table))) / length(union(colnames(seurat_table), colnames(hashtag_table)))

keep <- colnames(seurat_table) %in% colnames(hashtag_table)
seurat_table <- seurat_table[,keep]
keep <- colnames(hashtag_table) %in% colnames(seurat_table)
hashtag_table <- hashtag_table[,keep]

length(intersect(colnames(seurat_table), colnames(hashtag_table))) / length(union(colnames(seurat_table), colnames(hashtag_table)))

seurat_table[["Hash"]] <- CreateAssayObject( counts=hashtag_table )

rm(hashtag_table)
rm(keep)

seurat_table <- NormalizeData(seurat_table, assay = "Hash", normalization.method = "CLR")
seurat_table <- ScaleData(seurat_table, assay = "Hash")

#####
#filter data

seurat_table[["percent.mt"]] <- PercentageFeatureSet(seurat_table, pattern = "^mt-")

x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_hline( yintercept=100, color="red" ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="steelblue" ) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="pink" ) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=10, color="red" ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="orchid4" ) +
  theme_bw() +
  xlab("")

pdf( paste0("figures/qc_violin_prefilter_",date,".pdf" ), height=3, width=5 )
plot_grid( x,y,z, ncol=3 )
dev.off()

# VlnPlot(seurat_table, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_table <- subset(seurat_table, subset = percent.mt < 5)
seurat_table <- subset(seurat_table, subset = nFeature_RNA > 200)

x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_hline( yintercept=100, color="red" ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="steelblue" ) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="pink" ) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=10, color="red" ) +
  # geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="orchid4" ) +
  theme_bw() +
  xlab("")

pdf( paste0("figures/qc_violin_postfilter_",date,".pdf" ), height=3, width=5 )
plot_grid( x,y,z, ncol=3 )
dev.off()

#####
#remove doublets and assign hash identities

pdf( paste0("figures/qc_hash_",date,".pdf" ), height=5, width=8 )
for (i in 1:11) {
  x <- ggplot( as.data.frame(t(seurat_table@assays$Hash[i,])) ) +
    geom_histogram( aes( x=seurat_table@assays$Hash[i,] ),
                    fill=cmap[i], color="black" ) +
    geom_vline( xintercept=2 ) +
    geom_vline( xintercept=3 ) +
    geom_vline( xintercept=4 ) +
    # xlim(0,6) +
    theme_classic() +
    ggtitle( rownames(seurat_table@assays$Hash[i,]) ) +
    xlab( rownames(seurat_table@assays$Hash[i,]) )
  print(x)
}
dev.off()

png(paste0("figures/hash_pairs_plot_",date,".png"), height=1500, width=1500 )
pairs( t(seurat_table@assays$RNA[,]), cex=0.2, pch=20 )
dev.off()

sum( seurat_table@assays$Hash[1,] > 2 )
sum( seurat_table@assays$Hash[2,] > 2 )
sum( seurat_table@assays$Hash[3,] > 2 )
sum( seurat_table@assays$Hash[4,] > 1.5 )
sum( seurat_table@assays$Hash[5,] > 1.75 )
sum( seurat_table@assays$Hash[6,] > 2 )
sum( seurat_table@assays$Hash[7,] > 1.5 )
sum( seurat_table@assays$Hash[8,] > 1.5 )
sum( seurat_table@assays$Hash[9,] > 2 )
sum( seurat_table@assays$Hash[10,] > 2.5 )
sum( seurat_table@assays$Hash[11,] > 1.5 )

keep <- vector(mode="logical", length=ncol(seurat_table))
for( i in 1:ncol(seurat_table) ) {
  a <- seurat_table@assays$Hash[1,i] > 2
  b <- seurat_table@assays$Hash[2,i] > 2
  c <- seurat_table@assays$Hash[3,i] > 2
  d <- seurat_table@assays$Hash[4,i] > 1.5
  e <- seurat_table@assays$Hash[5,i] > 1.75
  f <- seurat_table@assays$Hash[6,i] > 2
  g <- seurat_table@assays$Hash[7,i] > 1.5
  h <- seurat_table@assays$Hash[8,i] > 1.5
  j <- seurat_table@assays$Hash[9,i] > 2
  k <- seurat_table@assays$Hash[10,i] > 2.5
  l <- seurat_table@assays$Hash[11,i] > 1.5
  if ( sum(a+b+c+d+e+f+g+h+j+k+l)!=1 ) {
    keep[i] = F
  } else {
    keep[i] = T
  }
}
seurat_table <- seurat_table[,keep]

#assign hash ident
ident <- vector(mode="logical", length=ncol(seurat_table))
for( i in 1:ncol(seurat_table) ) {
  if( seurat_table@assays$Hash[1,i] > 2 ) {
    ident[i] <- "wt_r1"
  } else if ( seurat_table@assays$Hash[2,i] > 2 ) {
    ident[i] <- "wt_r2"
  } else if ( seurat_table@assays$Hash[3,i] > 2 ) {
    ident[i] <- "wt_r3"
  } else if ( seurat_table@assays$Hash[4,i] > 1.5 ) {
    ident[i] <- "aire_r1"
  } else if ( seurat_table@assays$Hash[5,i] > 1.75 ) {
    ident[i] <- "aire_r2"
  } else if ( seurat_table@assays$Hash[6,i] > 2 ) {
    ident[i] <- "aire_r3"
  } else if ( seurat_table@assays$Hash[7,i] > 1.5 ) {
    ident[i] <- "spib_r1"
  } else if ( seurat_table@assays$Hash[8,i] > 1.5 ) {
    ident[i] <- "spib_r2"
  } else if ( seurat_table@assays$Hash[9,i] > 2 ) {
    ident[i] <- "sox8_r1"
  } else if ( seurat_table@assays$Hash[10,i] > 2.5 ) {
    ident[i] <- "sox8_r2"
  } else if ( seurat_table@assays$Hash[11,i] > 1.5 ) {
    ident[i] <- "wt_r4"
  } else {
    ident[i] <- "na"
  }
}

names(ident) <- colnames(seurat_table)
seurat_table@meta.data$hash.ident <- factor(ident)

ident <- vector(mode="logical", length=ncol(seurat_table))
for( i in 1:ncol(seurat_table) ) {
  if( seurat_table@assays$Hash[1,i] > 2 ) {
    ident[i] <- "wt"
  } else if ( seurat_table@assays$Hash[2,i] > 2 ) {
    ident[i] <- "wt"
  } else if ( seurat_table@assays$Hash[3,i] > 2 ) {
    ident[i] <- "wt"
  } else if ( seurat_table@assays$Hash[4,i] > 1.5 ) {
    ident[i] <- "aire"
  } else if ( seurat_table@assays$Hash[5,i] > 1.75 ) {
    ident[i] <- "aire"
  } else if ( seurat_table@assays$Hash[6,i] > 2 ) {
    ident[i] <- "aire"
  } else if ( seurat_table@assays$Hash[7,i] > 1.5 ) {
    ident[i] <- "spib"
  } else if ( seurat_table@assays$Hash[8,i] > 1.5 ) {
    ident[i] <- "spib"
  } else if ( seurat_table@assays$Hash[9,i] > 2 ) {
    ident[i] <- "sox8"
  } else if ( seurat_table@assays$Hash[10,i] > 2.5 ) {
    ident[i] <- "sox8"
  } else if ( seurat_table@assays$Hash[11,i] > 1.5 ) {
    ident[i] <- "wt"
  } else {
    ident[i] <- "na"
  }
}

names(ident) <- colnames(seurat_table)
seurat_table@meta.data$geno.ident <- factor(ident, levels=c("wt","aire","spib","sox8"))


#####
#normalize data
seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)

nfeatures = 2000
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)

top10 <- head(VariableFeatures(seurat_table), 10)

plot1 <- VariableFeaturePlot(seurat_table)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf( paste0("figures/variable_feature_postfilter_",date,".pdf" ), height=4, width=6 )
plot1
plot2
dev.off()

#####
#run dimensionality reduction
all.genes <- rownames(seurat_table)
seurat_table <- ScaleData(seurat_table, features = VariableFeatures(object = seurat_table))
seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))

sink( "text_outputs/top_PC_genes_1-30_postfilter.txt", append=F )
print(seurat_table[["pca"]], dims = 1:30, nfeatures = 5)
sink()

pdf( paste0("figures/PC_dim_loading_postfilter_",date,".pdf" ), height=40, width=12 )
VizDimLoadings(seurat_table, dims = 1:30, reduction = "pca")
dev.off()

DimPlot(seurat_table, reduction = "pca")

pdf( paste0("figures/PC_dim_heatmap_postfilter_",date,".pdf" ), height=40, width=12 )
DimHeatmap(seurat_table, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

seurat_table <- JackStraw(seurat_table, num.replicate = 100, dims=50)
seurat_table <- ScoreJackStraw(seurat_table, dims = 1:50)

pdf( paste0("figures/PC_jackstraw_",date,".pdf" ), height=5, width=12 )
JackStrawPlot(seurat_table, dims = 1:50)
dev.off()

pdf( paste0("figures/PC_elbow_postfilter_",date,".pdf" ), height=5, width=6 )
ElbowPlot(seurat_table, ndims=50)
dev.off()

#####
# cluster and umap
seurat_table <- FindNeighbors(seurat_table, dims = 1:30)
seurat_table <- FindClusters(seurat_table, resolution = 1.2)
seurat_table <- RunUMAP(seurat_table, dims = 1:30, seed.use = 1)

pdf( paste0("figures/umap_postfilter_",date,".pdf" ), height=5, width=7 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, cols=cmap)
dev.off()

#####
#remove contaminating cells
keep <- seurat_table@active.ident!=19 & seurat_table@active.ident!=19 #t,endothelial 
seurat_table <- seurat_table[,keep]

#####
#combine homologous clusters
keep <- seurat_table@active.ident==5 | seurat_table@active.ident==9 | seurat_table@active.ident==14 | seurat_table@active.ident==10
seurat_table <- seurat_table[,keep]

#####
#plot per hash umap
keep1 <- seurat_table@meta.data$hash.ident %in% c( "wt_r1","wt_r2","wt_r3", "wt_r4" )
keep2 <- seurat_table@meta.data$hash.ident %in% c( "aire_r1","aire_r2", "aire_r3" )
keep3 <- seurat_table@meta.data$hash.ident %in% c( "spib_r1","spib_r2" )
keep4 <- seurat_table@meta.data$hash.ident %in% c( "sox8_r1","sox8_r2" )
w <- DimPlot(seurat_table, group.by="geno.ident", pt.size=0.5, order=F, cols=cmap )
x <- DimPlot(seurat_table[,keep1], group.by="hash.ident", pt.size=0.5, order=F, cols=cmap[c(3,4,5,6)] )
y <- DimPlot(seurat_table[,keep2], group.by="hash.ident", pt.size=0.5, order=F, cols=cmap )
z <- DimPlot(seurat_table[,keep3], group.by="hash.ident", pt.size=0.5, order=F, cols=cmap )
a <- DimPlot(seurat_table[,keep4], group.by="hash.ident", pt.size=0.5, order=F, cols=cmap )

pdf( paste0("figures/umap_perGenotype_",date,".pdf" ), height=5, width=30 )
plot_grid(w,x,y,z,a, ncol=5)
dev.off()

w <- DimPlot(seurat_table, group.by="geno.ident", pt.size=0.5, order=F, cols=cmap )
x <- DimPlot(seurat_table[,keep1], group.by="geno.ident", pt.size=0.5, order=F, cols="gray" )
y <- DimPlot(seurat_table[,keep2], group.by="geno.ident", pt.size=0.5, order=F, cols="gray" )
z <- DimPlot(seurat_table[,keep3], group.by="geno.ident", pt.size=0.5, order=F, cols="gray" )
a <- DimPlot(seurat_table[,keep4], group.by="geno.ident", pt.size=0.5, order=F, cols="gray" )

pdf( paste0("figures/umap_perGenotype_monochrome_",date,".pdf" ), height=5, width=30 )
plot_grid(w,x,y,z,a, ncol=5)
dev.off()

w <- DimPlot(seurat_table, group.by="geno.ident", pt.size=0.5, order=F, cols=cmap )
x <- DimPlot(seurat_table[,keep1], pt.size=0.5, order=F, cols=cmap ) + theme(legend.position="none")
y <- DimPlot(seurat_table[,keep2], pt.size=0.5, order=F, cols=cmap ) + theme(legend.position="none")
z <- DimPlot(seurat_table[,keep3], pt.size=0.5, order=F, cols=cmap ) + theme(legend.position="none")
a <- DimPlot(seurat_table[,keep4], pt.size=0.5, order=F, cols=cmap ) + theme(legend.position="none")

pdf( paste0("figures/umap_perGenotype_perCluster_",date,".pdf" ), height=5, width=30 )
plot_grid(w,x,y,z,a, ncol=5)
dev.off()

#####
#calculate fractional representation of clusters in neonate vs adult

source("functions/calculate_cluster_fraction.R")

cluster_frac_list <- calculate_cluster_fraction(n_hash=11, 
                                                names_hash=unique(seurat_table@meta.data$hash.ident),
                                                n_cluster=length(unique(seurat_table@active.ident)),
                                                names_cluster=unique(seurat_table@active.ident),
                                                n_geno=4,
                                                names_geno=c("wt","aire","spib","sox8"),
                                                hash_path=seurat_table@meta.data$hash.ident,
                                                cluster_path=seurat_table@active.ident)

pdf( paste0("figures/cluster_abundance_boxplot_",date,".pdf"), height=4, width=12 )
ggplot( cluster_frac_list[[1]], aes( x=cluster, y=frac, fill=geno ) ) +
  # geom_jitter( aes(color=genotype), position=position_dodge( width = 0.75 ) ) +
  geom_boxplot( width=0.4, position=position_dodge( width = 0.75 ), lwd=0.3, outlier.shape = NA ) +
  # stat_summary(geom = "pointrange", aes(color=genotype), shape=3,
  #              fun.data = mean_se, position=position_dodge( width = 0.75 ) ) +
  scale_fill_manual( values=cmap ) +
  scale_color_manual( values=cmap ) +
  geom_vline(xintercept = seq(1.5,17.5, by=1), lwd=0.25, lty=3 ) +
  theme_bw() +
  theme( axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1), legend.position=c(0.9,0.6) ) +
  xlab("") +
  ylab( "Fraction of total Pdpn-CD104- MEClo")
dev.off()

pdf( paste0("figures/cluster_abundance_barplot_",date,".pdf"), height=12, width=6 )
ggplot(cluster_frac_list[[2]], aes(fill=cluster, y=frac, x=genotype, color="black")) + 
  geom_bar(position="fill", stat="identity", width=0.65) +
  scale_fill_manual( values=cmap ) +
  scale_color_manual( values="black") +
  xlab("Cluster") +
  ylab("Fraction of cells") +
  theme_few() +
  labs(color="") +
  theme(legend.position="right")
dev.off()

#####
#differential density plot

h=c(10,10)

densityWT <- kde2d( x=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("wt_r1","wt_r2","wt_r3","wt_r4"),"UMAP_1"], 
                    y=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("wt_r1","wt_r2","wt_r3","wt_r4"),"UMAP_2"],
                    h=h,
                    n=100 )
densityAire <- kde2d( x=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("aire_r1","aire_r2","aire_r3"),"UMAP_1"], 
                      y=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("aire_r1","aire_r2","aire_r3"),"UMAP_2"],
                      h=h,
                      n=100 )

densityAll <- kde2d( x=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("wt_r1","wt_r2","wt_r3","wt_r4","aire_r1","aire_r2","aire_r3"),"UMAP_1"], 
                     y=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("wt_r1","wt_r2","wt_r3","wt_r4","aire_r1","aire_r2","aire_r3"),"UMAP_2"],
                     h=h,
                     n=100 )

densityDiff <- densityWT$z - densityAire$z

xy = xy.coords( x=seurat_table@reductions$umap@cell.embeddings[,"UMAP_1"],
                y=seurat_table@reductions$umap@cell.embeddings[,"UMAP_2"] )
select = is.finite(xy$x) & is.finite(xy$y)
x = cbind(xy$x, xy$y)[select, ]
mkBreaks = function(u) u - diff(range(u))/(length(u) - 1)/2
xbin = cut(x[, 1], mkBreaks(densityAll$x), labels = FALSE)
ybin = cut(x[, 2], mkBreaks(densityAll$y), labels = FALSE)
dens = densityDiff[cbind(xbin, ybin)]
dens[is.na(dens)] = 0

pdf( paste0("figures/wt_aire_diff_density_",date,".pdf"), height=5, width=7 )
ggplot( as.data.frame( seurat_table@reductions$umap@cell.embeddings ),
        aes( x=seurat_table@reductions$umap@cell.embeddings[,"UMAP_1"],
             y=seurat_table@reductions$umap@cell.embeddings[,"UMAP_2"] ) ) +
  geom_point( aes( color=oob_squish(dens, range=c(-0.003,0.003)) ), size=0.5 ) +
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
  scale_color_viridis_c() +
  labs( color="Density differential\n(WT - Aire)" ) +
  xlab("UMAP-1") +
  ylab("UMAP-2")
dev.off()

#####
#save work
# saveRDS(seurat_table, file = "Rdata/meclo2_seurat_table_2021-11-23.rds")

library(tidyverse)
library(Seurat)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(BuenColors)
library(ggrepel)

setwd( "C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_100821" )
seurat_table <- readRDS( "Rdata/meclo2_seurat_table_2021-11-23.rds" )

date="20211124"

cmap2 = jdb_color_maps
cmap2 <- as.vector(cmap2)
cmap2 <- c(cmap2,"darkseagreen4")
names(cmap2) <- c("Aire-expressing","Tuft","Microfold","Enterocyte/hepatocyte","Neuroendocrine","Immature","Ionocyte","Aire-deficient","Ptf1a+ ductal",
                  "Keratinocyte","TA","Basal","Ciliated","Lung","perinatal cTEC","Muscle","adult cTEC","Tuft2")
cmap2["adult cTEC"] <- "darkorchid4"
cmap2["Lung"] <- "#993299"
cmap2["Neuroendocrine"] <- "#79adba"
cmap2["Immature"] <- "#f9838b"
cmap <- cmap2

#####
#viz data

feature <- "Chat"

a <- DimPlot(seurat_table, reduction="umap", label=T, cols=cmap, repel=T ) + theme(legend.position="none")
x <- FeaturePlot(seurat_table, features = c(feature), pt.size=1, order=T ) + theme(legend.position="none")
y <- VlnPlot( seurat_table, features=c(feature), cols=cmap ) + theme(legend.position="none")

# pdf( paste0("figures/umap_",feature, "_",date,".pdf"), height=4, width=5 )
plot_grid(a,x,y, ncol=3)
# dev.off()

#####
#analyze cluster composition
cluster.markers <- FindMarkers(seurat_table, ident.1="Aire-deficient", min.pct=0.1, logfc.threshold=1, only.pos=F, test.use="bimod")
head(cluster.markers, n = 20)

#####
#analyze cluster composition

idx1 <- (seurat_table@meta.data$hash.ident %in% c("wt_r1","wt_r2","wt_r3","wt_r4") & seurat_table@active.ident %in% c("Microfold") )
idx2 <- (seurat_table@meta.data$hash.ident %in% c("spib_r1","spib_r2") & seurat_table@active.ident %in% c("Microfold") )
idx3 <- (seurat_table@meta.data$hash.ident %in% c("sox8_r1","sox8_r2") & seurat_table@active.ident %in% c("Microfold") )
idx4 <- (seurat_table@meta.data$hash.ident %in% c("aire_r1","aire_r2","aire_r3") & seurat_table@active.ident %in% c("Microfold") )

set.seed(12345)

cluster.markers.spib <- FindMarkers(seurat_table@assays$RNA, 
                               cells.1=colnames(seurat_table)[idx1], 
                               cells.2=colnames(seurat_table)[idx2], 
                               only.pos=F, logfc.threshold=0, min.pct=0, test.use="bimod")
cluster.markers.sox8 <- FindMarkers(seurat_table@assays$RNA, 
                                cells.1=colnames(seurat_table)[idx1], 
                                cells.2=colnames(seurat_table)[idx3], 
                                only.pos=F, logfc.threshold=0, min.pct=0, test.use="bimod")
cluster.markers.aire <- FindMarkers(seurat_table@assays$RNA, 
                                cells.1=colnames(seurat_table)[idx1], 
                                cells.2=colnames(seurat_table)[idx4], 
                                only.pos=F, logfc.threshold=0.0, min.pct=0, test.use="bimod")

cluster.markers <- cluster.markers.spib
signature <- as.vector(read.delim(paste0("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures/20210924_Mcell_sig_logfc1_pct10_qval0.01.txt"), header=F, sep="\t")$V1)
color_idx <- ifelse(rownames(cluster.markers)%in%signature,"1","0")
label_idx <- ifelse( abs(cluster.markers$avg_log2FC) > 1, rownames(cluster.markers), "" )

ggplot( cluster.markers, aes(x=-avg_log2FC, y=-log10(p_val_adj) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_point( data=cluster.markers[signature,], aes( x=-cluster.markers[signature,"avg_log2FC"], 
                                                     y=-log10(cluster.markers[signature,"p_val_adj"]) ), 
              color="purple" ) +
  # geom_text_repel(label=label_idx ) +
  scale_color_manual( values=c("gray","purple") ) +
  geom_hline( yintercept=2, lty="dashed" ) +
  geom_vline( xintercept=1, lty="dashed" ) +
  geom_vline( xintercept=-1, lty="dashed" ) +
  theme_few() +
  xlim(-4,4) +
  ylim(0,25) +
  theme(legend.position="none") +
  xlab( paste0("log2 FC") ) +
  ylab( paste0("-log10 p-value") ) +
  ggtitle("")

#####
#make specific volcanos

cluster.markers <- FindMarkers(seurat_table, ident.1="Aire-deficient", min.pct=0, logfc.threshold=0.2, only.pos=F, test.use="bimod")
head(cluster.markers, n = 20)

label_idx <- ifelse(rownames(cluster.markers)%in%c("Irf8","H2-Aa","H2-Ab1","Cd74","H2-Oa","H2-Ob","H2-Eb1","H2-DMb2","H2-Eb2"), rownames(cluster.markers),"")
color_idx <- ifelse( cluster.markers$avg_log2FC > 1 & -log10(cluster.markers$p_val_adj) > 2, "1", "0"  )
color_idx <- ifelse( cluster.markers$avg_log2FC < -1 & -log10(cluster.markers$p_val_adj) > 2, "2", color_idx  )

pdf( paste0("figures/aire-deficient_volcano_",date,".pdf"), height=4, width=5 )
ggplot( cluster.markers, aes(x=avg_log2FC, y=-log10(p_val_adj) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_text_repel( label=label_idx, max.overlaps=50, force=20 ) +
  geom_hline( yintercept=2, lty="dashed" ) +
  geom_vline( xintercept=1, lty="dashed" ) +
  geom_vline( xintercept=-1, lty="dashed" ) +
  scale_color_manual( values=c("gray","navyblue","firebrick") ) +
  xlim(-4,4) +
  theme_few() +
  xlab( paste0("log2 FC, Aire-deficient vs all") ) +
  ylab( paste0("-log10 p-value") ) +
  theme( legend.position="none" )
dev.off()

cluster.markers <- cluster.markers.aire
signature <- as.vector(read.delim(paste0("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures/20210924_Skin, keratinized_sig_logfc1_pct10_qval0.01.txt"), header=F, sep="\t")$V1)
color_idx <- ifelse(rownames(cluster.markers)%in%signature,"1","0")
label_idx <- ifelse( -(cluster.markers$avg_log2FC) > 0 & rownames(cluster.markers)%in%signature, "up", "" )
label_idx <- ifelse( -(cluster.markers$avg_log2FC) < 0 & rownames(cluster.markers)%in%signature, "down", label_idx )

pdf( paste0("figures/keratinocyte_volcano_",date,".pdf"), height=4, width=5 )
ggplot( cluster.markers, aes(x=-avg_log2FC, y=-log10(p_val) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_point( data=cluster.markers[signature,], aes( x=-cluster.markers[signature,"avg_log2FC"], 
                                                     y=-log10(cluster.markers[signature,"p_val"]) ), 
              color="purple" ) +
  # geom_text_repel(label=label_idx ) +
  scale_color_manual( values=c("gray","purple") ) +
  geom_hline( yintercept=2, lty="dashed" ) +
  geom_vline( xintercept=1, lty="dashed" ) +
  geom_vline( xintercept=-1, lty="dashed" ) +
  theme_few() +
  scale_x_continuous( limits=c(-4,4), oob=squish ) +
  scale_y_continuous( limits=c(0,20), oob=squish ) +
  theme(legend.position="none") +
  xlab( paste0("log2 FC, Aire+/+ vs Aire-/- keratinocyte mTEC") ) +
  ylab( paste0("-log10 p-value") ) +
  ggtitle("keratinocyte") +
  annotate( geom="text", label=sum(label_idx=="up"), x=3.5, y=18, color="purple" ) +
  annotate( geom="text", label=sum(label_idx=="down"), x=-3.5, y=18, color="purple" )
dev.off()


#####
#build cluster tree
seurat_table <- BuildClusterTree(seurat_table, dims=T)
PlotClusterTree(seurat_table)

#####
#rename clusters
seurat_table <- RenameIdents( object=seurat_table,
                              "0"="Tuft",
                              "1"="Tuft",
                              "2"="Basal",
                              "3"="Tuft",
                              "4"="Aire-expressing",
                              "5"="Tuft",
                              "6"="Aire-deficient",
                              "7"="Tuft",
                              "8"="Neuroendocrine",
                              "9"="Keratinocyte",
                              "10"="Immature",
                              "11"="TA",
                              "12"="Tuft",
                              "13"="Lung",
                              "14"="Aire-expressing",
                              "15"="Microfold",
                              "16"="Aire-deficient",
                              "17"="Enterocyte/hepatocyte",
                              "18"="Ciliated",
                              "19"="Muscle")

seurat_table@active.ident <- factor(seurat_table@active.ident, 
                                    levels=c( "Tuft","Aire-expressing",
                                              "Enterocyte/hepatocyte","Aire-deficient",
                                              "Neuroendocrine","Keratinocyte","Immature",
                                              "TA","Microfold","Ciliated","Muscle","Lung","Basal" ) )

#####
#analyze cluster distribution in spib/sox8 and aire
source("functions/calculate_cluster_fraction.R")

cluster_frac_list <- calculate_cluster_fraction(n_hash=11, 
                           names_hash=unique(seurat_table@meta.data$hash.ident),
                           n_cluster=length(unique(seurat_table@active.ident)),
                           names_cluster=unique(seurat_table@active.ident),
                           n_geno=4,
                           names_geno=c("wt","aire","spib","sox8"),
                           hash_path=seurat_table@meta.data$hash.ident,
                           cluster_path=seurat_table@active.ident)

temp <- cluster_frac_list[[1]]
temp <- temp[temp$geno%in%c("wt","spib","sox8")==T & temp$cluster%in%c("Microfold"),]

pdf( paste0("figures/microfold_cluster_abundance_boxplot_",date,".pdf"), height=4, width=3 )
ggplot( temp, aes( x=cluster, y=frac, fill=geno ) ) +
  geom_boxplot( width=0.4, position=position_dodge( width = 0.75 ), lwd=0.3, outlier.shape = NA ) +
  geom_jitter( position=position_dodge( width = 0.75 ) ) +
  # stat_summary(geom = "pointrange", aes(color=genotype), shape=3,
  #              fun.data = mean_se, position=position_dodge( width = 0.75 ) ) +
  scale_fill_manual( values=c("gray","steelblue","royalblue") ) +
  scale_color_manual( values=c("gray","steelblue","royalblue") ) +
  ylim(0,max(temp$frac)) +
  # geom_vline(xintercept = seq(1.5,17.5, by=1), lwd=0.25, lty=3 ) +
  theme_bw() +
  theme( axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="right" ) +
  xlab("") +
  ylab( "Fraction of total Pdpn-CD104- MEClo")
dev.off()

temp <- cluster_frac_list[[2]]
temp <- temp[temp$geno%in%c("wt","spib","sox8")==T & temp$cluster%in%c("Aire-expressing","TA","Immature")==F,]

pdf( paste0("figures/microfold_cluster_abundance_barplot_",date,".pdf"), height=12, width=6 )
ggplot(temp, aes(fill=cluster, y=frac, x=genotype, color="black")) + 
  geom_bar(position="fill", stat="identity", width=0.65) +
  scale_fill_manual( values=cmap ) +
  scale_color_manual( values="black") +
  xlab("Cluster") +
  ylab("Fraction of cells") +
  theme_few() +
  labs(color="") +
  theme(legend.position="right")
dev.off()

temp <- cluster_frac_list[[1]]
temp <- temp[temp$geno%in%c("wt","aire")==T,]

pdf( paste0("figures/aire_cluster_abundance_boxplot_",date,".pdf"), height=4, width=12 )
ggplot( temp, aes( x=cluster, y=frac, fill=geno ) ) +
  geom_boxplot( width=0.4, position=position_dodge( width = 0.75 ), lwd=0.3, outlier.shape = NA ) +
  geom_jitter( position=position_dodge( width = 0.75 ) ) +
  # stat_summary(geom = "pointrange", aes(color=genotype), shape=3,
  #              fun.data = mean_se, position=position_dodge( width = 0.75 ) ) +
  scale_fill_manual( values=c("darkgreen","purple3") ) +
  scale_color_manual( values=c("darkgreen","purple3") ) +
  ylim(0,max(temp$frac)) +
  geom_vline(xintercept = seq(1.5,17.5, by=1), lwd=0.25, lty=3 ) +
  theme_few() +
  theme( axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="right" ) +
  xlab("") +
  ylab( "Fraction of total Pdpn-CD104- MEClo")
dev.off()

temp <- cluster_frac_list[[2]]
temp <- temp[temp$geno%in%c("wt","aire")==T & temp$cluster%in%c("Aire-expressing","TA","Immature")==F,]

pdf( paste0("figures/aire_cluster_abundance_barplot_",date,".pdf"), height=12, width=5 )
ggplot(temp, aes(fill=cluster, y=frac, x=genotype, color="black")) + 
  geom_bar(position="fill", stat="identity", width=0.65) +
  scale_fill_manual( values=cmap ) +
  scale_color_manual( values="black") +
  xlab("Cluster") +
  ylab("Fraction of cells") +
  theme_few() +
  labs(color="") +
  theme(legend.position="right")
dev.off()

temp <- cluster_frac_list[[1]]
temp <- temp[temp$geno%in%c("wt","aire","spib","sox8")==T & temp$cluster%in%c("Aire-expressing","TA","Immature")==F,]
order_cluster <- c("Microfold","Aire-deficient","Enterocyte/hepatocyte","Neuroendocrine","Keratinocyte","Ciliated","Muscle","Lung","Basal","Tuft")

pdf( paste0("figures/all_geno_cluster_abundance_barplot_perHash_",date,".pdf"), height=6, width=7 )
ggplot(temp, aes(fill=factor(cluster, levels=order_cluster), y=frac, x=factor(hash, levels=c("wt_r2","wt_r1","wt_r3","wt_r4",
                                                                                             "spib_r1","spib_r2",
                                                                                             "sox8_r1","sox8_r2",
                                                                                             "aire_r2","aire_r1","aire_r3")), color="black")) + 
  geom_bar(position="fill", stat="identity", width=0.65) +
  scale_fill_manual( values=cmap ) +
  scale_color_manual( values="black") +
  geom_vline(xintercept=c(4.5,6.5,8.5)) +
  xlab("Cluster") +
  ylab("Fraction of cells") +
  theme_few() +
  labs(color="") +
  theme(legend.position="right")
dev.off()

temp <- cluster_frac_list[[1]]

