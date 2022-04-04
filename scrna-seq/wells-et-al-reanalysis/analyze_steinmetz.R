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
library(future)

plan("multiprocess", workers = 4)

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/steinmetz_scrna")

date="2021-09-16"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value

#####
#read data

seurat_table <- read.delim("data/GSM4085087_aireTrace.csv", header=T, sep=",")

seurat_table <- CreateSeuratObject(seurat_table, min.cells = 1, project="AireLinTrace")

#####
#normalize data and find variable features
seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)

nfeatures = 1000
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)

top10 <- head(VariableFeatures(seurat_table), 10)

plot1 <- VariableFeaturePlot(seurat_table)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf( paste0("figures/variable_feature_cleaned_",date,".pdf" ), height=4, width=6 )
plot1
plot2
dev.off()

#####
#run dimensionality reduction

variable.genes <- VariableFeatures(object = seurat_table)
seurat_table <- ScaleData(seurat_table, features = variable.genes)
seurat_table <- RunPCA(seurat_table, features = variable.genes)

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

pdf( paste0("figures/PC_jackstraw_cleaned_",date,".pdf" ), height=5, width=12 )
JackStrawPlot(seurat_table, dims = 1:50)
dev.off()

pdf( paste0("figures/PC_elbow_",date,".pdf" ), height=5, width=6 )
ElbowPlot(seurat_table, ndims=50)
dev.off()

#####
# cluster and umap

seurat_table <- FindNeighbors(seurat_table, dims = 1:20)
seurat_table <- FindClusters(seurat_table, resolution = 2)

head(Idents(seurat_table), 100)

seurat_table <- RunUMAP(seurat_table, dims = 1:20, seed.use = 123)

pdf( paste0("figures/umap_labeled_",date,".pdf" ), height=5, width=7 )
DimPlot(seurat_table, reduction="umap", label=T, repel=T )
dev.off()

#####
#remove contaminating cells
keep <- seurat_table@active.ident!=15 & seurat_table@active.ident!=16 #t,myeloid 
seurat_table <- seurat_table[,keep]

#####
#combine homologous clusters
keep <- seurat_table@active.ident==5 | seurat_table@active.ident==9 | seurat_table@active.ident==14 | seurat_table@active.ident==10
seurat_table <- seurat_table[,keep]

#####
#save work
# saveRDS(seurat_table, file = "Rdata/steinmetz_seurat_table_2021-09-16.rds")

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

setwd( "C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/steinmetz_scrna" )
seurat_table <- readRDS( "Rdata/steinmetz_seurat_table_2021-09-16.rds" )

date="2021-09-17"
cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value

#####
#viz data

FeaturePlot(seurat_table, features = c("GFP"), order=T ) + theme(legend.position="none")

a <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Grhl1"), pt.size=1 ) + theme(legend.position="none")
b <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Foxa2"), pt.size=1 ) + theme(legend.position="none")
c <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Hnf4g"), pt.size=1 ) + theme(legend.position="none")
d <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Spib"), pt.size=1 ) + theme(legend.position="none")
e <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Sox8"), pt.size=1 ) + theme(legend.position="none")
f <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Foxj1"), pt.size=1 ) + theme(legend.position="none")
g <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Myog"), pt.size=1 ) + theme(legend.position="none")
h <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Pou2f3"), pt.size=1 ) + theme(legend.position="none")
i <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Aire"), pt.size=1 ) + theme(legend.position="none")
j <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Krt5"), pt.size=1 ) + theme(legend.position="none")
k <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Ccl21a"), pt.size=1 ) + theme(legend.position="none")
l <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Pdpn"), pt.size=1 ) + theme(legend.position="none")

pdf( paste0("figures/zsGreen_vs_postAireTF_correlation_grid_2021-09-16.pdf"), height=12, width=16 )
plot_grid(a,b,c,d,e,f,g,h,i,j,k,l)
dev.off()


DimPlot(seurat_table, reduction="umap", label=T, repel=T )
VlnPlot( joint_table, features=c("Hnf4a"), cols=cmap ) + theme(legend.position="none")

DimPlot(seurat_table, reduction="umap", label=T, repel=T )
cluster.markers <- FindMarkers(seurat_table, ident.1=c("25"), min.pct=0.1, logfc.threshold=1, only.pos=T)
head(cluster.markers, n = 20)

feature <- "Cd24"

a <- DimPlot(seurat_table, reduction="umap", label=F, cols=cmap, repel=T ) + theme(legend.position="none")
x <- FeaturePlot(seurat_table, features = c(feature), pt.size=1, order=T ) + theme(legend.position="none")
y <- VlnPlot( seurat_table, features=c(feature), cols=cmap ) + theme(legend.position="none")

# pdf( paste0("figures/2020-10-12/umap/umap_",feature, "_2020-10-12.pdf"), height=4, width=12 )
plot_grid(a,x,y, ncol=3)
# dev.off()

#####
#analyze post-aire sig corr

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures")
idx <- list.files()[grepl("20210924_", list.files())]
sig_list <- list()
for (i in idx ) {
  sig <- read.delim(i,header=T, sep="\t")
  sig <- sig[,1]
  sig_name <- substr(i, 10, nchar(i)-30)
  sig_list[[sig_name]] <- sig
}
setwd( "C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/steinmetz_scrna" )

for (sig_name in names(sig_list)){
  seurat_table <- AddModuleScore( seurat_table, features=as.data.frame(sig_list[[sig_name]]), name=sig_name )
}

pdf( paste0("figures/postAireSig_umap_grid_2021-10-10a_newsig.pdf"), height=12, width=16 )
FeaturePlot(seurat_table, features=c(names(seurat_table@meta.data)[6:23]), order=T, min.cutoff="q50", max.cutoff="q100", cols=c("lightgray","#3F007D"))
dev.off()

a <- FeaturePlot(seurat_table,order=T, features=c("ciliated1") ) + theme(legend.position="none")
b <- FeaturePlot(seurat_table,order=T, features=c("goblet1") ) + theme(legend.position="none")
c <- FeaturePlot(seurat_table,order=T, features=c("ionocyte1") ) + theme(legend.position="none")
d <- FeaturePlot(seurat_table,order=T, features=c("lung_basal1") ) + theme(legend.position="none")
e <- FeaturePlot(seurat_table,order=T, features=c("mcell1") ) + theme(legend.position="none")
f <- FeaturePlot(seurat_table,order=T, features=c("muscle1") ) + theme(legend.position="none")
g <- FeaturePlot(seurat_table,order=T, features=c("neuroendocrine1") ) + theme(legend.position="none")
h <- FeaturePlot(seurat_table,order=T, features=c("ptf1a_ductal1") ) + theme(legend.position="none")
i <- FeaturePlot(seurat_table,order=T, features=c("skin_basal1") ) + theme(legend.position="none")
j <- FeaturePlot(seurat_table,order=T, features=c("skin_keratinized1") ) + theme(legend.position="none")
k <- FeaturePlot(seurat_table,order=T, features=c("Tuft11") ) + theme(legend.position="none")
l <- FeaturePlot(seurat_table,order=T, features=c("Tuft21") ) + theme(legend.position="none")
m <- FeaturePlot(seurat_table,order=T, features=c("gutliver1") ) + theme(legend.position="none")
n <- FeaturePlot(seurat_table,order=T, features=c("Aire") ) + theme(legend.position="none")
o <- FeaturePlot(seurat_table,order=T, features=c("Krt10") ) + theme(legend.position="none")
p <- FeaturePlot(seurat_table,order=T, features=c("GFP") ) + theme(legend.position="none")

pdf( paste0("figures/postAireSig_umap_grid_2021-09-16.pdf"), height=16, width=16 )
plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p, ncol=4)
dev.off()

a <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Ciliated1"), pt.size=1 ) + theme(legend.position="none")
b <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Goblet1"), pt.size=1 ) + theme(legend.position="none")
c <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Gut1"), pt.size=1 ) + theme(legend.position="none")
d <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Ionocyte1"), pt.size=1 ) + theme(legend.position="none")
e <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Lung..basal1"), pt.size=1 ) + theme(legend.position="none")
f <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Mcell1"), pt.size=1 ) + theme(legend.position="none")
g <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Muscle1"), pt.size=1 ) + theme(legend.position="none")
h <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Neuroendocrine1"), pt.size=1 ) + theme(legend.position="none")
i <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Ptf1a..ductal1"), pt.size=1 ) + theme(legend.position="none")
j <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Skin..basal1"), pt.size=1 ) + theme(legend.position="none")
k <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Skin..keratinized1"), pt.size=1 ) + theme(legend.position="none")
l <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Tuft11"), pt.size=1 ) + theme(legend.position="none")
m <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Tuft21"), pt.size=1 ) + theme(legend.position="none")
n <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Aire"), pt.size=1 ) + theme(legend.position="none")
o <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Krt10"), pt.size=1 ) + theme(legend.position="none")
p <- FeatureScatter(seurat_table, feature1= c("GFP"), feature2=("Pdpn"), pt.size=1 ) + theme(legend.position="none")


pdf( paste0("figures/zsGreen_vs_postAireSig_correlation_grid_2021-09-16.pdf"), height=16, width=16 )
plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p, ncol=4)
dev.off()

pdf( paste0("figures/dimplot_2021-09-16.pdf"), height=6, width=7 )
DimPlot(seurat_table)
dev.off()

#####
#summarize sig/gfp relation

idx <- names(seurat_table@meta.data)[c(8:10,12:16,18:20,22:23)]

cutoff_table <- data.frame( cluster=idx, cutoffs=c(0.25,0.5,0.5,0.3,0.8,0.5,2,0.5,2,0.5,0.2,0.5,0.75) )
postaire_table <- matrix(nrow=length(idx), ncol=2)
colnames(postaire_table) <- c("aire_lin_neg","aire_lin_pos")
rownames(postaire_table) <- idx

for (i in idx ) {
  postaire_table[i,"aire_lin_pos"] <- sum( seurat_table@meta.data[,i] > cutoff_table[cutoff_table$cluster==i,2] & seurat_table@assays$RNA[,]["GFP",] > 1.6 )
  postaire_table[i,"aire_lin_neg"] <- sum( seurat_table@meta.data[,i] > cutoff_table[cutoff_table$cluster==i,2] & seurat_table@assays$RNA[,]["GFP",] < 1.6 )
}

features <- c("Aire","Krt10","Ccl21a")
cutoff_table <- data.frame( cluster=features, cutoffs=c(1,2,4) )
gene_table <- matrix(nrow=length(features), ncol=2)
colnames(gene_table) <- c("aire_lin_neg","aire_lin_pos")
rownames(gene_table) <- features

features <- c("Aire","Trp63","Grhl1","Foxa2","Hnf4g","Spib","Sox8","Foxj1","Myog","Pou2f3","Foxi1")
cutoff_table <- data.frame( cluster=features, cutoffs=c(1,1,1,1,1,2,1,2,2,1,1.5) )
gene_table <- matrix(nrow=length(features), ncol=2)
colnames(gene_table) <- c("aire_lin_neg","aire_lin_pos")
rownames(gene_table) <- features

for (i in features ) {
  gene_table[i,"aire_lin_pos"] <- sum( seurat_table@assays$RNA[,][i,] > cutoff_table[cutoff_table$cluster==i,2] & seurat_table@assays$RNA[,]["GFP",] > 1.6 )
  gene_table[i,"aire_lin_neg"] <- sum( seurat_table@assays$RNA[,][i,] > cutoff_table[cutoff_table$cluster==i,2] & seurat_table@assays$RNA[,]["GFP",] < 1.6 )
}

output_table <- as.data.frame(cbind(c(rownames(postaire_table),rownames(gene_table)),rbind(postaire_table, gene_table)))
output_table[,4] <- as.numeric(output_table[,2])+as.numeric(output_table[,3])
write_delim(output_table,"text_outputs/postaireSig_table_2021-09-17.txt",col_names=T,delim="\t")

output_table <- as.data.frame(cbind(rownames(gene_table),gene_table))
output_table[,4] <- as.numeric(output_table[,2])+as.numeric(output_table[,3])
write_delim(output_table,"text_outputs/postaireTF_table_2021-09-17.txt",col_names=T,delim="\t")

plot_table <- data.frame( lineage=c(rep("aire_lin_neg",times=length(idx)+3),rep("aire_lin_pos",times=length(idx)+3)),
                          cluster=factor(c(idx,features,idx,features), levels=c("Aire","Krt10","Ccl21a",idx)),
                          number=c(postaire_table[,1],gene_table[,1],postaire_table[,2],gene_table[,2]) )

plot_table <- data.frame( lineage=c(rep("aire_lin_neg",times=length(features)),rep("aire_lin_pos",times=length(features))),
                          cluster=factor(c(features,features), levels= c("Aire","Trp63","Foxj1","Hnf4g","Foxi1","Foxi2",
                                                                         "Foxa1","Foxa2","Foxa3","Spib","Sox8","Myog","Grhl1","Pou2f3")),
                          number=c(gene_table[,1],gene_table[,2]) )

pdf( "figures/cluster_postAire_TF_fraction_2021-10-10.pdf", height=5, width=14 )
ggplot( plot_table, aes(x=cluster, y=number, fill=lineage) ) +
  geom_bar(position="fill",stat="identity", color="black", width=0.6) +
  geom_text(aes(label=number), position='fill') +
  theme_bw() +
  xlab("") +
  ylab("Fraction of cells in cluster") +
  scale_fill_manual(values=c("darkgreen","maroon"))
dev.off()

#####
#analyze cluster composition
cluster.markers <- FindMarkers(seurat_table, ident.1=c(3), min.pct=0.01, logfc.threshold=0.5, only.pos=T)
head(cluster.markers, n = 20)
# write_delim( as.data.frame(rownames(cluster.markers)), "scratch/test.txt", delim="\t", col_names=F )

#####
#build cluster tree
seurat_table <- BuildClusterTree(seurat_table, dims=T)
PlotClusterTree(seurat_table)

#####
#assign cluster names
seurat_table <- RenameIdents( object=seurat_table,
                              "33"="Ciliated",
                              "30"="Ciliated",
                              "29"="Lung2",
                              "31"="Neuro",
                              "18"="Skin",
                              "25"="Lung1",
                              "15"="Gut/liver",
                              "22"="Mcell",
                              "3"="Aire-stage",
                              "4"="Aire-stage",
                              "17"="Transitional",
                              "32"="Aire-adjacent",
                              "24"="Tuft",
                              "26"="Tuft",
                              "19"="TA",
                              "28"="TA",
                              "23"="TA",
                              "6"="cTEC1",
                              "7"="cTEC1",
                              "14"="cTEC1",
                              "2"="cTEC1",
                              "34"="cTEC1",
                              "1"="cTEC1",
                              "9"="cTEC1",
                              "8"="cTEC1",
                              "11"="cTEC1",
                              "5"="cTEC2",
                              "12"="cTEC2",
                              "13"="cTEC2",
                              "27"="cTEC2",
                              "0"="cTEC2",
                              "10"="cTEC2",
                              "16"="cTEC2",
                              "21"="cTEC2",
                              "20"="cTEC2")

FeaturePlot(seurat_table, features = c("Spib"), pt.size=1, order=T ) + theme(legend.position="none")
DimPlot(seurat_table, pt.size=1, order=T, label=T, cols=cmap )

seurat_table@active.ident <- factor(seurat_table@active.ident, 
                                    levels=c( "Ciliated","Lung2","Neuro","Skin","Lung1",
                                              "Gut/liver","Mcell","Aire-stage","Transitional",
                                              "Aire-adjacent","Tuft","TA","cTEC1","cTEC2") )
