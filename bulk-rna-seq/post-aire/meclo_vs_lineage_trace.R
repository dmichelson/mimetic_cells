#####
#read in lineage trace data

library(tidyverse)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/2020-2-8_MEC_lineage_tracing")

date = "2021-09-24"

#####
#read in data
big_table <- read.delim( "data/GSM3146499_RFPn_hi_1.txt", header=T, sep=" " )
rownames(big_table) <- big_table[,1]
big_table <- select(big_table,2)

for (i in dir(path="data/.") ) {
  x <- read.delim( paste("data/",i,sep=""), header=T, sep=" " )
  rownames(x) <- x[,1]
  x <- select(x,2)
  big_table <- data.frame( big_table, x )
  rm(x)
}

head(big_table)

big_table <- select(big_table, -c(RFPn_hi_1.1))

gene_list <- as.data.frame(read_delim("gene_list_biomart_2-11-2020.txt", delim="\t", col_names=T))

remove <- gene_list[,2]==""
gene_list_culled <- gene_list[!remove,]
remove <- gene_list_culled[,1]==""
gene_list_culled <- gene_list_culled[!remove,]
gene_list_culled <- gene_list_culled[!duplicated(gene_list_culled[,1]),]
gene_list_culled <- gene_list_culled[!duplicated(gene_list_culled[,2]),]
gene_list_culled <- na.omit(gene_list_culled)

keep <- rownames(big_table)%in%gene_list_culled[,1]
big_table_cleaned <- big_table[keep,]
rownames(big_table_cleaned) <- gene_list_culled[,2]

#####
library(edgeR)

group <- factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4))
y <- DGEList( counts=big_table_cleaned, group=group, genes=rownames(big_table_cleaned) )

keep <- filterByExpr(y, min.count=10)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y,design)

qlf.RFPnHi_vs_RFPnLo <- glmQLFTest(fit,coef=2)
qlf.RFPnHi_vs_RFPpHi <- glmQLFTest(fit,coef=3)
qlf.RFPnHi_vs_RFPpLo <- glmQLFTest(fit,coef=4)
qlf.RFPnLo_vs_RFPpLo <- glmQLFTest(fit,contrast=c(0,1,0,-1))
qlf.RFPpHi_vs_RFPpLo <- glmQLFTest(fit,contrast=c(0,0,1,-1))
qlf.RFPnLo_vs_RFPpHi <- glmQLFTest(fit,contrast=c(0,1,-1,0))
qlf.RFPpLo_vs_RFPnLo <- glmQLFTest(fit,contrast=c(0,-1,0,1))

topTags(qlf.RFPnLo_vs_RFPpLo, n=50)

#####
setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2020-08-21_mcell_meclo_rnaseq/Batch096_Michelson_Mouse")

norm_table <- read.delim( "gene_expression_1_adjust.gct", sep="\t", header=T, skip=2, row.names=1 )
norm_table <- norm_table[,c(5:7,14:19)]
count_table <- read.delim( "genes_count_table_QC_passed.csv", sep=",", header=T, row.names=1 )
count_table <- count_table[,c(4:6,13:18)]

sample_names <- colnames(norm_table)

setwd("meclo")

#####
# de analysis
i = "PDPNp.CD104p.MEClo.Thy"
j = "PDPNn.CD104n.MEClo.Thy"

idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == i
ctrl_idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == j

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

qlf.PDPNnCD104nLO_vsPDPNpCD104pLO <- qlf

#####

qlf.1 <- qlf.RFPpLo_vs_RFPnLo$table
qlf.2 <- qlf.PDPNnCD104nLO_vsPDPNpCD104pLO$table

ident.1 <- "Aire+ vs Aire-"
ident.2 <- "Pdpn-CD104- vs Pdpn+CD104+"

keep <- rownames(qlf.1) %in% rownames(qlf.2)
qlf.1 <- qlf.1[keep,]
keep <- rownames(qlf.2) %in% rownames(qlf.1)
qlf.2 <- qlf.2[keep,]

qlf.1 <- qlf.1[order(rownames(qlf.1)),]
qlf.2 <- qlf.2[order(rownames(qlf.2)),]

library(Hmisc)
temp <- rcorr(x=qlf.1$logFC,y=qlf.2$logFC)

pdf( paste0("figures/meclo_marker_lineage_FCFC_rcorr_",date,".pdf"), height=5, width=6 )
label_idx <- ifelse( qlf.1$logFC < -3 & qlf.2$logFC < -3, rownames(qlf.1), "" )
ggplot( qlf.1, aes( x=qlf.1$logFC, y=qlf.2$logFC ) ) +
  geom_hline( yintercept = 0 ) +
  geom_vline( xintercept = 0 ) +
  geom_point( aes() ) +
  # geom_smooth( method="loess" ) +
  # geom_text_repel( label=label_idx ) +
  theme_bw() +
  xlab(ident.1) +
  ylab(ident.2) +
  theme(legend.position="none") +
  scale_x_continuous( limits=c(-6,6), oob=squish ) +
  scale_y_continuous( limits=c(-6,8), oob=squish ) +
  scale_color_continuous( low=c("white"), high=c("black") ) +
  annotate( x=5, y=7, geom="text", label=paste0( "r=",substr(temp$r[1,2],1,4) ) )
dev.off()

#####
#signatures

path_name="C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures"
idx <- list.files(path=path_name)[grepl(".txt", list.files(path=path_name))]

pdf( paste0("figures/meclo_marker_lineage_FCFC_allSigs_",date,".pdf"), height=5, width=6 )
for (i in idx ) {
  sig <- read.delim(paste0(path_name,"/",i),header=T, sep="\t")
  sig <- sig[,1]
  sig <- sig[sig%in%rownames(qlf.1)]
  sig_name <- substr(i, 10, nchar(i)-30)

  color_idx <- ifelse( rownames(qlf.1) %in% sig, "1", "0" )

  p <- ggplot( qlf.1, aes( x=qlf.1$logFC, y=qlf.2$logFC ) ) +
    geom_hline( yintercept = 0 ) +
    geom_vline( xintercept = 0 ) +
    geom_point( aes( color=color_idx ) ) +
    geom_point( data=as.data.frame(qlf.1[sig,]), aes( x=qlf.1[sig,"logFC"], y=qlf.2[sig,"logFC"] ), color="purple" ) +
    theme_bw() +
    xlab(ident.1) +
    ylab(ident.2) +
    theme(legend.position="none") +
    scale_color_manual( values=c("gray","purple") ) +
    ggtitle(sig_name)
  print(p)
}
dev.off()
