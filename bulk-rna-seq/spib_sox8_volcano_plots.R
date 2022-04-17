library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(pcaMethods)
library(ggrepel)
library(scales)

date = "2022-02-06"

#####
#sox8
setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2021-11-02_sox8_dnmeclo_rnaseq/Batch119_Dan_Mouse")

norm_table <- read.delim( "Batch119_Mouse_Dan_1_adjust.gct", sep="\t", header=T, skip=2, row.names=1 )
norm_table <- norm_table[,2:7]
count_table <- read.delim( "Genes_count_table.tsv", sep="\t", header=T, row.names=1 )
count_table <- count_table[,1:6]

sample_names <- colnames(norm_table)

i = "Sox8WT"
j = "Sox8KO"

idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 15 ) == i
ctrl_idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 15 ) == j

group <- factor(c(rep(1,times=sum(idx)),rep(2,times=sum(ctrl_idx))))
y <- DGEList( counts=cbind( count_table[,idx], count_table[,ctrl_idx]), 
              group=group, 
              genes=rownames(count_table) )

keep <- filterByExpr(y, min.count=5)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y,design)

qlf <- glmQLFTest(fit,coef=2)

qlf$table <- cbind( qlf$table, p.adjust( qlf$table$PValue, method="BH" ) )
colnames(qlf$table)[5] = "q_val"

qlf.sox8 <- qlf

#####
#spib

setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2021-08-19_spib_ckm_rnaseq/Batch116_Michelson_Mouse")

norm_table <- read.delim( "Batch116_Mouse_Michelson_1_adjust.gct", sep="\t", header=T, skip=2, row.names=1 )
norm_table <- norm_table[,2:13]
count_table <- read.delim( "Genes_count_table.tsv", sep="\t", header=T, row.names=1 )
count_table <- count_table[,1:12]

sample_names <- colnames(norm_table)

i = "SpibWT.DN.MEClo"
j = "SpibKO.DN.MEClo"

idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 3 ) == i
ctrl_idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 3 ) == j

# idx <- substr( colnames(count_table), 1, 6 ) == i
# ctrl_idx <- substr( colnames(count_table), 1, 6 ) == j

group <- factor(c(rep(1,times=sum(idx)),rep(2,times=sum(ctrl_idx))))
y <- DGEList( counts=cbind( count_table[,idx], count_table[,ctrl_idx]), 
              group=group, 
              genes=rownames(count_table) )

keep <- filterByExpr(y, min.count=5)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y,design)

qlf <- glmQLFTest(fit,coef=2)

qlf$table <- cbind( qlf$table, p.adjust( qlf$table$PValue, method="BH" ) )
colnames(qlf$table)[5] = "q_val"

qlf.spib <- qlf

#####
#plot
setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2021-11-02_sox8_dnmeclo_rnaseq/Batch119_Dan_Mouse")
qlf <- qlf.sox8

i = "SpibWT.DN.MEClo"
j = "SpibKO.DN.MEClo"

i = "Sox8WT"
j = "Sox8KO"

# mcell_sig <- read.delim("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_scATAC/2020-06-05_organ_sig/haber_nature_2017_mcellSig.txt", header=T, sep="\t")
# mcell_sig <- unique(c( mcell_sig[,1],mcell_sig[,2] ))

# mcell_sig <- read.delim("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-06-07_sig/2021-06-07_mcell_sig_logfc1_pct10_onlypos.txt", header=F, sep="\t")
# mcell_sig <- mcell_sig$V1

mcell_sig <- read.delim("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures/20210924_Mcell_sig_logfc1_pct10_qval0.01.txt", header=F, sep="\t")
mcell_sig <- mcell_sig$V1

signature <- mcell_sig

label_idx <- ifelse( (rownames(qlf$table) %in% signature & qlf$table$logFC < -1) & (qlf$table$q_val < 0.2), rownames(qlf$table), "" )
color_idx <- vector( length=nrow(qlf$table) )
for (k in 1:nrow(qlf$table)) {
  if (rownames(qlf$table)[k] %in% signature){
    color_idx[k] <- "mcell"
  } else {
    color_idx[k] <- "not_mcell"
  }
}

pdf(paste0("figures/",i,"_vs_", j,"_volcano_mcell_bespoke_sig_lab_", date ,".pdf"), height=3, width=4.5)
ggplot( qlf$table, aes( x=qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_point( data=qlf$table[signature,], aes( x=qlf$table[signature,"logFC"], 
                                               y=-log10(qlf$table[signature,"PValue"]) ),
              color="purple" ) +
  # geom_hline(yintercept=1.3, linetype="dashed") +
  # geom_vline(xintercept=-1, linetype="dashed") +
  # geom_vline(xintercept=1, linetype="dashed") +
  geom_text_repel( label=label_idx, color="black", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab(paste0("log2 fold change, ", j, " vs ", i)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("purple","gray") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  ggtitle( paste0(j, " vs ", i) ) +
  scale_x_continuous(limits=c(-4,4), oob=squish) +
  scale_y_continuous(limits=c(0,6), oob=squish) +
  annotate( geom="text", label=sum(color_idx=="mcell"&qlf$table$logFC>0), x=4, y=0.5, color="purple" ) +
  annotate( geom="text", label=sum(color_idx=="mcell"&qlf$table$logFC<0), x=-4, y=0.5, color="purple" )
dev.off()
