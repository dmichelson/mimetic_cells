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

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620")

date="2021-11-06"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value

#####
#read data

gene_table <- Read10X( data.dir="filtered_feature_bc_matrix" )
seurat_table <- CreateSeuratObject( counts = gene_table, project = "PostAireMEClo", min.cells=3, min.features = 0 )

#####
#integrate hashes
hashtag_table <- t(read.table(file = "hash/umi_count/mat.txt", sep = "\t", header = TRUE, row.names = 1))
colnames(hashtag_table) <- paste0( colnames(hashtag_table), "-1" )

length(intersect(colnames(seurat_table), colnames(hashtag_table))) / length(union(colnames(seurat_table), colnames(hashtag_table)))

keep <- colnames(seurat_table) %in% colnames(hashtag_table)
seurat_table <- seurat_table[,keep]

keep <- colnames(hashtag_table) %in% colnames(seurat_table)
hashtag_table <- hashtag_table[,keep]

length(intersect(colnames(seurat_table), colnames(hashtag_table))) / length(union(colnames(seurat_table), colnames(hashtag_table)))

seurat_table[["Hash"]] <- CreateAssayObject( counts=hashtag_table )
seurat_table <- NormalizeData(seurat_table, assay = "Hash", normalization.method = "CLR")
seurat_table <- ScaleData(seurat_table, assay = "Hash")

#####
#filter data

seurat_table[["percent.mt"]] <- PercentageFeatureSet(seurat_table, pattern = "^mt-")

x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  # geom_hline( yintercept=200, color="red" ) +
  geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="steelblue" ) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="pink" ) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=25, color="red" ) +
  geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="orchid4" ) +
  theme_bw() +
  xlab("")

pdf( paste0("figures/qc_violin_",date,".pdf" ), height=3, width=5 )
plot_grid( x,y,z, ncol=3 )
dev.off()

VlnPlot(seurat_table, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_table <- subset(seurat_table, subset = percent.mt < 25)

#####
#viz hashes

colors <- palette("Tableau 10")

x <- ggplot( as.data.frame(t(seurat_table@assays$Hash[1,])) ) +
  geom_histogram( aes( x=seurat_table@assays$Hash[1,] ),
                  fill=colors[1], color="black" ) +
  geom_vline( xintercept=1.1 ) +
  theme_classic() +
  xlab( rownames(seurat_table@assays$Hash[1,]) )
y <- ggplot( as.data.frame(t(seurat_table@assays$Hash[2,])) ) +
  geom_histogram( aes( x=seurat_table@assays$Hash[2,] ),
                  fill=colors[2], color="black" ) +
  geom_vline( xintercept=1.1 ) +
  theme_classic() +
  xlab( rownames(seurat_table@assays$Hash[2,]) )
z <- ggplot( as.data.frame(t(seurat_table@assays$Hash[3,])) ) +
  geom_histogram( aes( x=seurat_table@assays$Hash[3,] ),
                  fill=colors[3], color="black" ) +
  geom_vline( xintercept=1.1 ) +
  theme_classic() +
  xlab( rownames(seurat_table@assays$Hash[3,]) )
a <- ggplot( as.data.frame(t(seurat_table@assays$Hash[4,])) ) +
  geom_histogram( aes( x=seurat_table@assays$Hash[4,] ),
                  fill=colors[4], color="black" ) +
  geom_vline( xintercept=3 ) +
  theme_classic() +
  xlab( rownames(seurat_table@assays$Hash[4,]) )
b <- ggplot( as.data.frame(t(seurat_table@assays$Hash[5,])) ) +
  geom_histogram( aes( x=seurat_table@assays$Hash[5,] ),
                  fill=colors[5], color="black" ) +
  geom_vline( xintercept=3 ) +
  theme_classic() +
  xlab( rownames(seurat_table@assays$Hash[5,]) )
c <- ggplot( as.data.frame(t(seurat_table@assays$Hash[6,])) ) +
  geom_histogram( aes( x=seurat_table@assays$Hash[6,] ),
                  fill=colors[6], color="black" ) +
  geom_vline( xintercept=3 ) +
  theme_classic() +
  xlab( rownames(seurat_table@assays$Hash[6,]) )

pdf( paste0("figures/qc_hash_",date,".pdf" ), height=5, width=8 )
plot_grid( x,y,z,a,b,c, ncol=3)
dev.off()

sum( seurat_table@assays$Hash[1,] > 1.1 )
sum( seurat_table@assays$Hash[2,] > 1.1 )
sum( seurat_table@assays$Hash[3,] > 1.1 )
sum( seurat_table@assays$Hash[4,] > 3 )
sum( seurat_table@assays$Hash[5,] > 3 )
sum( seurat_table@assays$Hash[6,] > 3 )

#####
#remove doublets

pdf( paste0("figures/qc_hash_doublets_",date,".pdf" ), height=8, width=8 )
pairs( t(seurat_table@assays$Hash[,]), cex=0.25, pch=16 )
dev.off()

keep <- vector(mode="logical", length=ncol(seurat_table))
for( i in 1:ncol(seurat_table) ) {
  a <- seurat_table@assays$Hash[1,i] > 1.1
  b <- seurat_table@assays$Hash[2,i] > 1.1
  c <- seurat_table@assays$Hash[3,i] > 1.1
  d <- seurat_table@assays$Hash[4,i] > 3
  e <- seurat_table@assays$Hash[5,i] > 3
  f <- seurat_table@assays$Hash[6,i] > 3
  if ( sum(a+b+c+d+e+f)!=1 ) {
    keep[i] = F
  } else {
    keep[i] = T
  }
}
seurat_table <- seurat_table[,keep]

#####
#normalize data and find variable features
seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)

nfeatures = 2000
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)

top10 <- head(VariableFeatures(seurat_table), 10)

plot1 <- VariableFeaturePlot(seurat_table)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf( paste0("figures/variable_feature_",date,".pdf" ), height=4, width=6 )
plot1
plot2
dev.off()

#####
#run dimensionality reduction
all.genes <- rownames(seurat_table)
seurat_table <- ScaleData(seurat_table, features = all.genes)

seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))

sink( "text_outputs/top_PC_genes_1-30.txt" )
print(seurat_table[["pca"]], dims = 1:30, nfeatures = 5)
sink()

pdf( paste0("figures/PC_dim_loading_",date,".pdf" ), height=40, width=12 )
VizDimLoadings(seurat_table, dims = 1:30, reduction = "pca")
dev.off()

DimPlot(seurat_table, reduction = "pca")

pdf( paste0("figures/PC_dim_heatmap_",date,".pdf" ), height=40, width=12 )
DimHeatmap(seurat_table, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

seurat_table <- JackStraw(seurat_table, num.replicate = 100, dims=50)
seurat_table <- ScoreJackStraw(seurat_table, dims = 1:50)

pdf( paste0("figures/PC_jackstraw_",date,".pdf" ), height=5, width=12 )
JackStrawPlot(seurat_table, dims = 1:50)
dev.off()

pdf( paste0("figures/PC_elbow_",date,".pdf" ), height=5, width=6 )
ElbowPlot(seurat_table, ndims=50)
dev.off()

#####
# cluster and umap

seurat_table <- FindNeighbors(seurat_table, dims = 1:40)
seurat_table <- FindClusters(seurat_table, resolution = 1.8)

head(Idents(seurat_table), 100)

seurat_table <- RunUMAP(seurat_table, dims = 1:40, seed.use = 12345)

pdf( paste0("figures/clusterAg_labeled_umap_",date,".pdf" ), height=5, width=7 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T, cols=cmap )
dev.off()

#####
#remove contaminating cells
keep <- seurat_table@active.ident!=23 & seurat_table@active.ident!=24 & seurat_table@active.ident!=25 & seurat_table@active.ident!=26 
seurat_table <- seurat_table[,keep]

#####
#assign hash ident
ident <- vector(mode="logical", length=ncol(seurat_table))
for( i in 1:ncol(seurat_table) ) {
  if( seurat_table@assays$Hash[1,i] > 1.1 ) {
    ident[i] <- "adult1"
  } else if ( seurat_table@assays$Hash[2,i] > 1.1 ) {
    ident[i] <- "adult2"
  } else if ( seurat_table@assays$Hash[3,i] > 1.1 ) {
    ident[i] <- "adult3"
  } else if ( seurat_table@assays$Hash[4,i] > 3 ) {
    ident[i] <- "neonate1"
  } else if ( seurat_table@assays$Hash[5,i] > 3 ) {
    ident[i] <- "neonate2"
  } else if ( seurat_table@assays$Hash[6,i] > 3 ) {
    ident[i] <- "neonate3"
  } else {
    ident[i] <- "na"
  }
}

names(ident) <- colnames(seurat_table)
seurat_table@meta.data$hash.ident <- factor(ident)
keep1 <- seurat_table@meta.data$hash.ident %in% c( "adult1","adult2","adult3" )
keep2 <- seurat_table@meta.data$hash.ident %in% c( "neonate2","neonate3" )
x <- DimPlot(seurat_table[,keep1], group.by="hash.ident", pt.size=1, order=T, cols=palette("Tableau 10") )
y <- DimPlot(seurat_table[,keep2], group.by="hash.ident", pt.size=1, order=T, cols=palette("Tableau 10") )

pdf( paste0("figures/neonate_vs_adult_umap_",date,".pdf" ), height=5, width=20 )
plot_grid(a,x,y, ncol=3)
dev.off()

#####
#calculate fractional representation of clusters in neonate vs adult

cluster_frac <- matrix( nrow=length(unique(seurat_table@active.ident)), ncol=5 )
colnames(cluster_frac) <- c("adult1", "adult2", "adult3","neonate2","neonate3")
rownames(cluster_frac) <- levels(seurat_table@active.ident)

for (i in levels(seurat_table@active.ident)) {
  print( paste0( "adult, cluster ", i, ":" ) )
  print( sum( keep1==T & seurat_table@active.ident==i ) / sum( keep1==T ) )
  cluster_frac[i,1] <- sum( seurat_table@meta.data$hash.ident=="adult1" & seurat_table@active.ident==i ) / sum( seurat_table@meta.data$hash.ident=="adult1" )
  cluster_frac[i,2] <- sum( seurat_table@meta.data$hash.ident=="adult2" & seurat_table@active.ident==i ) / sum( seurat_table@meta.data$hash.ident=="adult2" )
  cluster_frac[i,3] <- sum( seurat_table@meta.data$hash.ident=="adult3" & seurat_table@active.ident==i ) / sum( seurat_table@meta.data$hash.ident=="adult3" )
  # cluster_frac[i,4] <- sum( seurat_table@meta.data$hash.ident=="neonate1" & seurat_table@active.ident==i ) / sum( seurat_table@meta.data$hash.ident=="neonate1" )
  cluster_frac[i,4] <- sum( seurat_table@meta.data$hash.ident=="neonate2" & seurat_table@active.ident==i ) / sum( seurat_table@meta.data$hash.ident=="neonate2" )
  cluster_frac[i,5] <- sum( seurat_table@meta.data$hash.ident=="neonate3" & seurat_table@active.ident==i ) / sum( seurat_table@meta.data$hash.ident=="neonate3" )
  print( paste0( "neonate, cluster ", i, ":" ) )
  print( sum( keep2==T & seurat_table@active.ident==i ) / sum( keep2==T ) )
}

cluster_frac <- c( cluster_frac[,1], cluster_frac[,2], cluster_frac[,3], cluster_frac[,4], cluster_frac[,5] )
cluster_frac <- data.frame( cluster_frac, 
                            c(rep("adult", times=18*3), rep("neonate", times=18*2)), 
                            names(cluster_frac) )
colnames(cluster_frac) <- c("frac","age","cluster")

cluster_frac$cluster <- factor( cluster_frac$cluster, 
                                levels=(levels(seurat_table@active.ident)) )

pdf( "figures/cluster_abundance.pdf", height=4, width=8 )
ggplot( cluster_frac, aes( x=cluster, y=frac, fill=age ) ) +
  # geom_jitter( aes(color=age), position=position_dodge( width = 0.75 ) ) +
  geom_boxplot( width=0.4, position=position_dodge( width = 0.75 ), lwd=0.3, outlier.shape = NA ) +
  # stat_summary(geom = "pointrange", aes(color=age), shape=3,
  #              fun.data = mean_se, position=position_dodge( width = 0.75 ) ) +
  scale_fill_manual( values=cmap[c(1,3)] ) +
  scale_color_manual( values=cmap[c(1,3)] ) +
  geom_vline(xintercept = seq(1.5,17.5, by=1), lwd=0.25, lty=3 ) +
  theme_bw() +
  theme( axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1), legend.position=c(0.9,0.8) ) +
  xlab("") +
  ylab( "Fraction of total Pdpn-CD104- MEClo")
dev.off()

#####
#differential density plot

h=c(10,10)

densityNeonate <- kde2d( x=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("neonate2","neonate3"),"UMAP_1"], 
                         y=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("neonate2","neonate3"),"UMAP_2"],
                         h=h,
                         n=100 )
densityAdult <- kde2d( x=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("adult1","adult2","adult3"),"UMAP_1"], 
                       y=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("adult1","adult2","adult3"),"UMAP_2"],
                       h=h,
                       n=100 )

densityAll <- kde2d( x=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("neonate2","neonate3","adult1","adult2","adult3"),"UMAP_1"], 
                     y=seurat_table@reductions$umap@cell.embeddings[seurat_table@meta.data$hash.ident %in% c("neonate2","neonate3","adult1","adult2","adult3"),"UMAP_2"],
                     h=h,
                     n=100 )

densityDiff <- densityAdult$z - densityNeonate$z

xy = xy.coords( x=seurat_table@reductions$umap@cell.embeddings[,"UMAP_1"],
                y=seurat_table@reductions$umap@cell.embeddings[,"UMAP_2"] )
select = is.finite(xy$x) & is.finite(xy$y)
x = cbind(xy$x, xy$y)[select, ]
mkBreaks = function(u) u - diff(range(u))/(length(u) - 1)/2
xbin = cut(x[, 1], mkBreaks(densityAll$x), labels = FALSE)
ybin = cut(x[, 2], mkBreaks(densityAll$y), labels = FALSE)
dens = densityDiff[cbind(xbin, ybin)]
dens[is.na(dens)] = 0

pdf( "figures/adult_neonate_diff_density.pdf", height=5, width=7 )
ggplot( as.data.frame( seurat_table@reductions$umap@cell.embeddings ),
        aes( x=seurat_table@reductions$umap@cell.embeddings[,"UMAP_1"],
             y=seurat_table@reductions$umap@cell.embeddings[,"UMAP_2"] ) ) +
  geom_point( aes( color=oob_squish(dens, range=c(-0.003,0.003)) ), size=0.5 ) +
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
  scale_color_gradient2( low=cmap[3], mid="navajowhite2", high=cmap[1] ) +
  labs( color="Density differential\n(Adult - Neonate)" ) +
  xlab("UMAP-1") +
  ylab("UMAP-2")
dev.off()

#####
#combine homologous clusters
keep <- seurat_table@active.ident==2 | seurat_table@active.ident==4 | seurat_table@active.ident==5 | seurat_table@active.ident==7 | seurat_table@active.ident==8 | seurat_table@active.ident==10
seurat_table@active.ident[keep] <- 2

keep <- seurat_table@active.ident==0 | seurat_table@active.ident==3 | seurat_table@active.ident==16
seurat_table@active.ident[keep] <- 3

keep <- seurat_table@active.ident==1 | seurat_table@active.ident==13
seurat_table@active.ident[keep] <- 1

#####
#save work
# saveRDS(seurat_table, file = "Rdata/meclo_seurat_table_2021-09-16.rds")

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

setwd( "C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620" )
seurat_table <- readRDS( "Rdata/meclo_seurat_table_2021-09-16.rds" )

date="2021-11-06"

cmap2 = jdb_color_maps
cmap2 <- as.vector(cmap2)
cmap2 <- c(cmap2,"darkseagreen4")
names(cmap2) <- c("Aire-stage","Tuft1","Mcell","Gut/Liver","Neuroendocrine","Immature MEC","Ionocyte","Goblet","Ptf1a+ ductal",
                  "Skin, keratinized","TA MEC","Skin, basal","Ciliated","Lung","perinatal cTEC","Muscle","adult cTEC","Tuft2")
cmap2["adult cTEC"] <- "darkorchid4"
cmap2["Lung, basal"] <- "#993299"
cmap2["Neuroendocrine"] <- "#79adba"
cmap2["Immature MEC"] <- "#f9838b"
cmap <- cmap2

#####
#viz data

feature <- "Spib"

a <- DimPlot(seurat_table, reduction="umap", label=F, cols=cmap, repel=T ) + theme(legend.position="none")
x <- FeaturePlot(seurat_table, features = c(feature), pt.size=1, order=T ) + theme(legend.position="none")
y <- VlnPlot( seurat_table, features=c(feature), cols=cmap ) + theme(legend.position="none")

# pdf( paste0("figures/umap_",feature, ".pdf"), height=4, width=5 )
plot_grid(a,x,y, ncol=3)
# dev.off()

#####
#highlight parts of umap

feature <- "Dynlrb2"
# idx <- order(seurat_table@assays$RNA[feature,], decreasing=F)
pdf( paste0("figures/umap_",feature, "_",date,".pdf"), height=3, width=4 )
FeaturePlot( seurat_table, features = c(feature), pt.size=1, min.cutoff=0, max.cutoff=2, order=T, cols=c("lightgray","#3F007D"))
# VlnPlot( seurat_table, features=c(feature), cols=cmap, pt.size=0 ) + theme(legend.position="none")
dev.off()

#####
#analyze cluster composition
cluster.markers <- FindMarkers(seurat_table, ident.1="Gut/Liver", min.pct=0.01, logfc.threshold=0.5, only.pos=T)
head(cluster.markers, n = 20)
# write_delim( as.data.frame(rownames(cluster.markers)), "scratch/test.txt", delim="\t", col_names=F )

#####
#build cluster tree
seurat_table <- BuildClusterTree(seurat_table, dims=T)
PlotClusterTree(seurat_table)

#####
# make heatmap

features <- c( "Psmb11","Prss16","Castor1","Cxcl12","Ctsl","Ccl25","Tbata","Ly75",
               "Mki67","Birc5","Top2a","Cenpf","Hmgb2","Stmn1","Ccna2","Cenpe",
               "Trp63","Ccl21a","Itga6","Krt5","Stat1","Skint7","Mir205hg","Ptpn18",
               "Aire","Fezf2","Cd74","H2-Aa","Srgn","Ubd","H2-DMb2","S100a14",
               "Pou2f3","Chat","Avil","Trpm5","Alox5","Alox5ap","Dclk1","Ptprc",
               "Sox8","Spib","Ccl9","Ccl6","Ccl20","Anxa5","Tnfrsf11b","Tnfaip2",
               "Saa3","Apoa4","Ces2e","Guca2a","Reg3g","Lypd8","Nos2","Muc13",
               "Foxj1","Meig1","Spag8","Dynlrb2","Elof1","Ccdc153","Dnah12","Cfap44",
               "Myog", "Actc1", "Mymx", "Ckm","Cdkn1c","Meg3","Myl1","Klhl41",
               "Foxi1","Foxi2","Cftr","Atp6v1b1","Slc12a2","Ascl3","Slc4a11","Atp6v1g3",
               "Ptf1a","Kirrel2","Prss2","Clps","Gjb1","Csf2rb2","Rhov","Plvap",
               "Foxa2","Chgb","Ret","Snap25","Clca3b","Cacna2d1","Dnajc12","Car8",
               "Cp","Aqp4","Muc5b","Cxcl17","Wfdc2","Gsto1","Ltf","Serpinb11",
               "Serpinb2","Slco1a5","Gabrp","Spink5","Rptn","Gsta4","Oit1","Gm2a",
               "Grhl1","Dmkn","Il1f5","Dmkn","Krt10","Krtdap","Gltp","Lypd3",
               "Grhl1","Dmkn","Il1f5","Ivl","Asprv1","Cst6","Sbsn","Calm4")

anno <- vector(mode="character")

for( i in levels(seurat_table@active.ident) ) {
  temp <- rep( i, times=8 )
  anno <- c(anno, temp)
}

anno <- anno[!(anno==c("Tuft2"))]
anno[anno %in% c("Tuft1")] = "Tuft"

anno <- anno[!(anno==c("perinatal cTEC"))]
anno[anno %in% c("adult cTEC")] = "cTEC"

features <- cbind(features, anno)

features <- distinct( as.data.frame(features), features, .keep_all=T )

features <- as.matrix(features)

cells = WhichCells( seurat_table, 
                    idents=levels(seurat_table@active.ident)[levels(seurat_table@active.ident)%in%c("perinatal cTEC","Tuft1")==F],
                    downsample=50, seed=10 )

cluster_idx <- factor( seurat_table@active.ident[cells], levels=rev(levels(seurat_table@active.ident)) )

breakslist = seq(1,3,by=0.01)

ann_colors <- cmap[1:19]
ann_colors <- list( "cluster_idx"=ann_colors[levels(seurat_table@active.ident)] )

pdf( "figures/cluster_heatmap_test.pdf", height=7, width=12 )
pheatmap( seurat_table@assays$RNA@scale.data[features[,1],cells][,rev(order(cluster_idx))],
          show_rownames=T,
          show_colnames=F,
          cluster_cols=F,
          cluster_rows=F,
          annotation_col=as.data.frame(cluster_idx),
          annotation_colors = ann_colors,
          color=colorRampPalette((brewer.pal(n=9, name="Purples")))(length(breakslist)),
          breaks=breakslist,
          annotation_legend = T )
dev.off()

#####
#run all cluster de

library(future)
library(ggrepel)
plan("multiprocess", workers = 6)
cluster_de_list <- list()

# for (i in unique(seurat_table@active.ident)) {
  # if (i =="Gut/Liver") {
    # i="Gut"
    # cluster_de_list[[i]] <- read.delim( paste0("text_outputs/sig_de/cluster_de_list_cluster_",i,".txt"), header=T, sep="\t" )
    # rownames(cluster_de_list[[i]]) <- cluster_de_list[[i]][,1]
    # cluster_de_list[[i]] <- cluster_de_list[[i]][,-c(1)]
    # } else {
    # cluster_de_list[[i]] <- read.delim( paste0("text_outputs/sig_de/cluster_de_list_cluster_",i,".txt"), header=T, sep="\t" )
    # rownames(cluster_de_list[[i]]) <- cluster_de_list[[i]][,1]
    # cluster_de_list[[i]] <- cluster_de_list[[i]][,-c(1)]
  # }
# }

pdf( paste0("figures/cluster_specific_volcano_",date,".pdf"), height=4, width=5 )
for (i in unique(seurat_table@active.ident)) {
  print(i)
  cluster.markers <- FindMarkers(seurat_table, ident.1=i, min.pct=0, logfc.threshold=0.2, only.pos=F, test.use="bimod")
  head(cluster.markers, n = 20)
  cluster_de_list[[i]] <- cluster.markers
  cluster.markers <- cluster_de_list[[i]]
  if (i=="Gut/Liver") {
    cluster.markers <- cluster_de_list[[i]]
    i="Gut"
  } else { cluster.markers <- cluster_de_list[[i]] }
  
  label_idx <- ifelse( rownames(cluster.markers) %in% rownames(cluster.markers[1:15,]) & 
                         grepl("Rik", rownames(cluster.markers))==F &
                         grepl("Gm", rownames(cluster.markers))==F, rownames(cluster.markers), "" )
  color_idx <- ifelse( cluster.markers$avg_log2FC > 1 & -log10(cluster.markers$p_val_adj) > 2, "1", "0"  )
  color_idx <- ifelse( cluster.markers$avg_log2FC < -1 & -log10(cluster.markers$p_val_adj) > 2, "2", color_idx  )
  
  p <- ggplot(cluster.markers, aes(x=cluster.markers$avg_log2FC, y=-log10(cluster.markers$p_val_adj) )) +
    geom_point( aes(color=color_idx) ) +
    theme_bw() +
    geom_text_repel( label=label_idx, max.overlaps=50, force=20 ) +
    geom_hline( yintercept=2, lty="dashed" ) +
    geom_vline( xintercept=1, lty="dashed" ) +
    geom_vline( xintercept=-1, lty="dashed" ) +
    scale_color_manual( values=c("gray","navyblue","firebrick") ) +
    scale_x_continuous(limits=c(-4,4), oob=squish) +
    # scale_y_continuous(limits=c(0,300), oob=squish) +
    xlab( paste0("log2 FC, ", i," vs all") ) +
    ylab( paste0("-log10 p-value") ) +
    theme( legend.position="none" ) +
    ggtitle( paste0("cluster ",i) )
  print(p)
}
dev.off()

for (i in unique(seurat_table@active.ident)) {
  if (i =="Gut/Liver") {
    i="Gut"
    write_delim( as.data.frame( cbind(rownames(cluster_de_list[["Gut/Liver"]]),cluster_de_list[["Gut/Liver"]])), paste0("text_outputs/sig_de/cluster_de_list_cluster_",i,".txt"), col_names=T, delim="\t" )
  } else {
    write_delim( as.data.frame( cbind(rownames(cluster_de_list[[i]]),cluster_de_list[[i]])), paste0("text_outputs/sig_de/cluster_de_list_cluster_",i,".txt"), col_names=T, delim="\t" )
  }
}
 
levels(seurat_table) <- rev(levels(seurat_table))
 
pdf( "figures/all_clusters_dotplot.pdf", height=6, width=14 )
DotPlot(seurat_table, features=(c("Foxn1","Trp63","Aire","Pou2f3",
                                     "Spib","Sox8","Hnf4g","Hnf4a",
                                     "Foxj1","Trp73","Myog","Foxi1","Foxi2",
                                     "Ptf1a","Foxa2","Foxa3","Foxa1",
                                     "Spdef","Grhl1","Gata3")),
        cols=c("lightgray","#3F007D"), dot.min = 0.05, dot.scale=8 )
dev.off()


#####
# make signatures

intestinal_sig <- read.delim( "C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_scATAC/2020-06-05_organ_sig/haber_nature_2017_dropletIntestinalSig.txt",
                              sep="\t", header=T )

lung_sig <- read.delim( "C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_scATAC/2020-06-05_organ_sig/montoro_nature_2018_dropletLungSig.txt",
            sep="\t", header=T )

mcell_sig <- read.delim( "C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_scATAC/2020-06-05_organ_sig/haber_nature_2017_mcellSig.txt",
                   sep="\t", header=T )

skin_sig <- read.delim( "C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_scATAC/2020-06-05_organ_sig/joost_cellsystems_2016_skinSig.txt",
                         sep="\t", header=T )

hepatocyte_sig <- read.delim( "C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_scATAC/2020-06-05_organ_sig/aizarani_nature_2019_hepatocyteSig.txt",
                        sep="\t", header=F )

es_sig <- read.delim( "C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_scATAC/2020-06-05_organ_sig/ben-porath_natgen_esSig_2020-06-18.txt",
                      sep="\t", header=F )

muscle_sig <- read.delim( "C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/gene_lists/muscle_geneset_M40254.txt",
                      sep="\t", header=F )

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

intestinal_sig <- as.matrix(intestinal_sig)
eec_sig <- as.matrix(eec_sig)
lung_sig <- as.matrix(lung_sig)
mcell_sig <- as.matrix(mcell_sig)
skin_sig <- as.data.frame(skin_sig)
hepatocyte_sig <- str_to_sentence( hepatocyte_sig$V1 )
es_sig <- str_to_sentence( es_sig$V1 )
muscle_sig <- convertHumanGeneList( muscle_sig$V1[3:252] )

for (i in 1:ncol(skin_sig)) {
  skin_sig[,i] <- str_extract(skin_sig[,i], "^[^\\s]+")
}
  
seurat_table <- AddModuleScore(seurat_table, as.data.frame(intestinal_sig[intestinal_sig[,1]!="",1]), name="Intestine.Enteroendocrine")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(intestinal_sig[intestinal_sig[,2]!="",2]), name="Intestine.EnterocyteImmatureDistal")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(intestinal_sig[intestinal_sig[,3]!="",3]), name="Intestine.EnterocyteImmatureProximal")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(intestinal_sig[intestinal_sig[,4]!="",4]), name="Intestine.EnterocyteMatureDistal")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(intestinal_sig[intestinal_sig[,5]!="",5]), name="Intestine.EnterocyteMatureProximal")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(intestinal_sig[intestinal_sig[,9]!="",9]), name="Intestine.Goblet")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(intestinal_sig[intestinal_sig[,10]!="",10]), name="Intestine.Paneth")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(intestinal_sig[intestinal_sig[,15]!="",15]), name="Intestine.Tuft")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(mcell_sig[mcell_sig[,1]!="",1]), name="Intestine.Mcell.InVivo")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(mcell_sig[mcell_sig[,2]!="",2]), name="Intestine.McellInVitro")

seurat_table <- AddModuleScore(seurat_table, as.data.frame(lung_sig[lung_sig[,1]!="",1]), name="Lung.Basal")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(lung_sig[lung_sig[,2]!="",2]), name="Lung.Ciliated")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(lung_sig[lung_sig[,3]!="",3]), name="Lung.Club")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(lung_sig[lung_sig[,4]!="",4]), name="Lung.Goblet")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(lung_sig[lung_sig[,5]!="",5]), name="Lung.Ionocyte")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(lung_sig[lung_sig[,6]!="",6]), name="Lung.Neuroendocrine")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(lung_sig[lung_sig[,7]!="",7]), name="Lung.Tuft")

seurat_table <- AddModuleScore(seurat_table, as.data.frame(na.omit(skin_sig[,1])), name="Skin.Basal")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(na.omit(skin_sig[,3])), name="Skin.Differentiated.I")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(na.omit(skin_sig[,5])), name="Skin.Differentiated.II")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(na.omit(skin_sig[,7])), name="Skin.Keratinized.I")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(na.omit(skin_sig[,9])), name="Skin.Keratinized.II")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(na.omit(skin_sig[,11])), name="Skin.HairFollicle.I")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(na.omit(skin_sig[,13])), name="Skin.HairFollicle.II")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(na.omit(skin_sig[,15])), name="Skin.HairFollicle.III")

seurat_table <- AddModuleScore(seurat_table, as.data.frame(hepatocyte_sig), name="Hepatocyte")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(es_sig), name="Embryonic")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(muscle_sig), name="Muscle")
seurat_table <- AddModuleScore(seurat_table, as.data.frame(c(hepatocyte_sig,intestinal_sig[intestinal_sig[,5]!="",5])), name="GutLiver")

pdf( "figures/umap_signatures_muscleOnly.pdf", height=3.75, width=4 )
FeaturePlot(seurat_table, features = c("Lung.Ionocyte1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.2) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Skin.Keratinized.I1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.1) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Lung.Ciliated1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.1) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Lung.Basal1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.6) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Skin.Basal1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Intestine.EnterocyteMatureProximal1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.1) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Hepatocyte1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.1) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Intestine.McellInVitro1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.3) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Intestine.Tuft1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.25) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Intestine.EnterocyteMatureDistal1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Lung.Goblet1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.1) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Lung.Goblet1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.1) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Muscle1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.05) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("Lung.Neuroendocrine1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.2) + theme(legend.position="none")
FeaturePlot(seurat_table, features = c("GutLiver1"), cols=c("lightgray","#3F007D"), pt.size=0.75, order=T, min.cutoff=0.1) + theme(legend.position="none")
dev.off()


#####
#plot signature

# pdf( "figures/postAire_modules_validated_sig.pdf" )
FeaturePlot(seurat_table, features = c("Intestine.Enteroendocrine1", "Intestine.EnterocyteMatureProximal1","Intestine.Goblet1",
                                       "Intestine.Tuft1","Intestine.Mcell.InVivo1"),
            pt.size=0.5, sort.cell=T, min.cutoff = 0.1 )

FeaturePlot(seurat_table, features = c("Hepatocyte1"),pt.size=0.5, sort.cell=T, min.cutoff = 0.1 )

FeaturePlot(seurat_table, features = c("Skin.Basal1", "Skin.Differentiated.II1","Skin.Keratinized.I1",
                                       "Skin.HairFollicle.II1","skin.Langerhans1"),
            pt.size=0.5, sort.cell=T, min.cutoff = 0.1 )

FeaturePlot(seurat_table, features = c("Lung.Ciliated1","Lung.Club1","Lung.Goblet1",
                                       "Lung.Neuroendocrine1","Lung.Tuft1"),
            pt.size=0.5, sort.cell=T, min.cutoff = 0.1 )
# dev.off()
