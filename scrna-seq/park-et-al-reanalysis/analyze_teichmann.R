#####
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(harmony)

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/teichmann_scrna")

date="2021-11-06"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value

#####
#read data
# Convert("data/HTA08.v01.A05.Science_human_epi.gzip.h5ad", dest="data/epi.h5seurat", overwrite=T)

seurat_table <- LoadH5Seurat("data/epi.h5seurat")

#####
#analyze teichmann

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures")
idx <- list.files()[grepl("20210924_", list.files())]
sig_list <- list()
for (i in idx ) {
  sig <- read.delim(i,header=T, sep="\t")
  sig <- sig[,1]
  sig_name <- substr(i, 10, nchar(i)-30)
  sig_list[[sig_name]] <- sig
}
setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/teichmann_scrna")

convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

for (sig_name in names(sig_list)){
  sig <- convertMouseGeneList(sig_list[[sig_name]])
  seurat_table <- AddModuleScore( seurat_table, features=as.data.frame(sig), name=sig_name )
}

DimPlot(seurat_table, group.by="cell.types")

pdf( paste0("figures/",date,"_sig_overlay.pdf"), height=12, width=16)
FeaturePlot(seurat_table, features=c(names(seurat_table@meta.data)[12:29]), order=T, min.cutoff="q75", cols=c("lightgray","#3F007D"))
dev.off()

