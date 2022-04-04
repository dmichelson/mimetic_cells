library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(pcaMethods)
library(ggrepel)
library(scales)

setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2021-10-05_mcell_aire_rnaseq/Batch118_Michelson_Mouse")

date = "20211122"

norm_table <- read.delim( "Batch118_Mouse_Michelson_1_adjust.gct", sep="\t", header=T, skip=2, row.names=1 )
norm_table <- norm_table[,2:13]
count_table <- read.delim( "Genes_count_table.tsv", sep="\t", header=T, row.names=1 )
count_table <- count_table[,1:12]

sample_names <- colnames(norm_table)

#####
##QC

#look at counts data
pdf( paste0("figures/counts_per_sample_",date,".pdf"), height=3, width=5 )
hist(log10(colSums(count_table)), breaks=10, xlim=range(0,10))
abline( v=3.3010, col="red" )
dev.off()

#look at normalized data
pdf( paste0("figures/genes_per_sample_",date,".pdf"), height=3, width=5 )
hist(colSums( norm_table > 1 ), breaks=20, xlim=range(0,40000))
dev.off()

pdf( paste0("figures/detection_per_gene_",date,".pdf"), height=3, width=5 )
hist(rowSums( norm_table > 20 ), breaks=20)
abline(v=2, col="red")
dev.off()

#remove genes with normalized expression < threshold
idx <- rowSums( norm_table > 20 ) < 2
norm_table <- norm_table[!idx,]

# calculate basic stats for samples and genes
means <- rowMeans(norm_table)
vars <- apply(norm_table,1,var)
cvs <- apply(norm_table,1,sd)/rowMeans(norm_table)

# plot means vs cvs 
pdf( paste0("figures/cv_vs_mean_",date,".pdf"), height=5, width=5 )
plot(log(means),log(cvs))
dev.off()

#####
##initial clustering

#look at sample clustering
sample_cor_matrix <- cor( norm_table, method="pearson" )
hclust_samples <- hclust( as.dist( 1-sample_cor_matrix ), method="ward.D2" )

scaled_norm_table <- log2( norm_table )
scaled_norm_table <- scaled_norm_table - rowMeans(scaled_norm_table)

annotation_names = as.data.frame(sample_names)
rownames( annotation_names ) = rownames( sample_cor_matrix )

breaksList <- seq(0.2,1,by=0.01)

pdf( paste0("figures/samples_distance_heatmap_unclustered_", date, ".pdf"), height=5, width=6)
pheatmap( sample_cor_matrix, 
          cluster_rows=F,
          cluster_cols=F,
          show_annotation_row=T,
          labels_row = as.vector(sample_names), 
          show_rownames=T, 
          show_colnames=F, 
          main="Samples distance matrix", 
          annotation_legend = F, 
          scale="none",
          breaks=breaksList, 
          color=colorRampPalette((brewer.pal(name="Reds",n=9)))(length(breaksList)),
          border_color = NA
)
dev.off()

pdf( paste0("figures/samples_distance_heatmap_clustered_", date, ".pdf"), height=10, width=12)
pheatmap( sample_cor_matrix, 
          cluster_rows=hclust_samples,
          cluster_cols=hclust_samples,
          show_annotation_row=T,
          labels_row = as.vector(sample_names), 
          show_rownames=T, 
          show_colnames=F, 
          main="Samples distance matrix", 
          annotation_legend = F, 
          scale="none",
          breaks=breaksList, 
          color=colorRampPalette((brewer.pal(name="Reds",n=9)))(length(breaksList)),
          border_color = NA
)
dev.off()

#####
##pca analysis

#run PCA
pca_object <- pca( t(norm_table), method="svd", nPcs=10, scale="uv", center=T )

#examine elbow plot
pdf( paste0("figures/elbow_plot_", date, ".pdf" ), width=3.5, height=4 )
plot( pca_object@R2, xlab="PC#", ylab="R^2", main="Elbow plot" )
lines( pca_object@R2 )
dev.off()

#examine PC scores
pdf( paste0("figures/PCA_plot_", date, ".pdf" ), height=3.5, width=4 )
for (i in 1:9) {
  temp_plot <- ggplot( as.data.frame(pca_object@scores), aes( x=pca_object@scores[,i],y=pca_object@scores[,i+1] ) ) +
    geom_point( aes( color=sample_names ), size=4 ) +
    geom_text_repel( label=sample_names ) +
    theme_classic() +
    theme( legend.position="none" ) +
    xlab( paste0("PC",i, " (",pca_object@R2[i]*100,"%)") ) +
    ylab( paste0("PC",i+1, " (",pca_object@R2[i+1]*100,"%)") ) +
    theme( panel.border = element_rect(colour = "black", fill=NA, size=1) )
  print( temp_plot )
}
dev.off()

sample_names[1:12] <- c("Microfold.WT.1","Microfold.WT.2","Microfold.WT.3","Microfold.WT.4",
                       "Microfold.AireKO.1","Microfold.AireKO.2","Microfold.AireKO.3","Microfold.AireKO.4",
                       "mTEC.1","mTEC.2","mTEC.3","mTEC.4")

pdf( paste0("figures/PCA_plot_prettified_", date, ".pdf" ), height=3, width=4 )
i=1
temp_plot <- ggplot( as.data.frame(pca_object@scores), aes( x=pca_object@scores[,i],y=pca_object@scores[,i+1] ) ) +
  geom_point( aes( fill= factor( substr( sample_names, 1, nchar(sample_names)-2 ), 
                                 levels=c("Microfold.WT","mTEC","Microfold.AireKO") ) ), size=4, shape=21 ) +
  geom_text_repel( label=sample_names ) +
  theme_classic() +
  theme( legend.position="none" ) +
  xlab( paste0("PC",i, " (",pca_object@R2[i]*100,"%)") ) +
  ylab( paste0("PC",i+1, " (",pca_object@R2[i+1]*100,"%)") ) +
  scale_fill_manual( values=brewer.pal(n=4, name="Set1") ) +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1) )
print( temp_plot )
dev.off()


#find top loadings for PCs
pca_genes <- matrix( ncol=10, nrow=50 )
colnames( pca_genes ) = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
for (i in 1:10) {
  temp <- sort( pca_object@loadings[, paste0("PC",i) ], decreasing=T )
  pca_genes[,i] = names( temp[1:50] )
}

#write top loadings to file
write_delim( as.data.frame( pca_genes ), paste0( "de_analysis/top_pca_genes_per_component_", date, ".txt" ), delim="\t", col_names=T )

#find bottom loadings for PCs
pca_genes <- matrix( ncol=10, nrow=50 )
colnames( pca_genes ) = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
for (i in 1:10) {
  temp <- sort( pca_object@loadings[, paste0("PC",i) ], decreasing=F )
  pca_genes[,i] = names( temp[1:50] )
}

#write bottom loadings to file
write_delim( as.data.frame( pca_genes ), paste0( "de_analysis/bottom_pca_genes_per_component_", date, ".txt" ), delim="\t", col_names=T )

#####
##de analysis for thymus

# run edgeR pipeline
i = "Mcell.WT.Thy"
j = "Mcell.AireKO.Thy"

idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == i
ctrl_idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == j

group <- factor(c(rep(1,times=sum(idx)),rep(2,times=sum(ctrl_idx))))
y <- DGEList( counts=cbind( count_table[,idx], count_table[,ctrl_idx]), 
              group=group, 
              genes=rownames(count_table) )

keep <- filterByExpr(y, min.count=20)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y,design)

qlf <- glmQLFTest(fit,coef=2)

qlf$table <- cbind( qlf$table, p.adjust( qlf$table$PValue, method="BH" ) )
colnames(qlf$table)[5] = "q_val"

qlf.aireko <- qlf

write_delim( cbind( rownames(qlf$table),qlf$table ), paste0("de_analysis/",i,"_vs_", j,"_DE_table_", date ,".txt"), col_names=T, delim="\t" )

i = "Mcell.WT.Thy"
j = "GP2n.MEC.Thy"

idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == i
ctrl_idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == j

group <- factor(c(rep(1,times=sum(idx)),rep(2,times=sum(ctrl_idx))))
y <- DGEList( counts=cbind( count_table[,idx], count_table[,ctrl_idx]), 
              group=group, 
              genes=rownames(count_table) )

keep <- filterByExpr(y, min.count=20)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y,design)

qlf <- glmQLFTest(fit,coef=2)

qlf$table <- cbind( qlf$table, p.adjust( qlf$table$PValue, method="BH" ) )
colnames(qlf$table)[5] = "q_val"

qlf.mcell <- qlf

write_delim( cbind( rownames(qlf$table),qlf$table ), paste0("de_analysis/",i,"_vs_", j,"_DE_table_", date ,".txt"), col_names=T, delim="\t" )

#plot exp-exp
qlf <- qlf.aireko
i = "Mcell.WT.Thy"
j = "Mcell.AireKO.Thy"
idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == i
ctrl_idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == j

label_idx <- ifelse( rownames(norm_table) %in% rownames(topTags(qlf, n=50)) &
                       !grepl("Rik", rownames(norm_table)), rownames(norm_table), "" )
color_idx <- ifelse( rownames(norm_table) %in% rownames(topTags(qlf, n=50)), "1", "0" )

pdf(paste0("figures/",i,"_vs_", j,"_exp-exp_", date ,".pdf"), height=8, width=9)
ggplot( norm_table, aes( x=rowMeans(log2(norm_table[,ctrl_idx])), y=rowMeans(log2(norm_table[,idx])) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab( paste0( "log2 expression, ",j ) ) +
  ylab( paste0( "log2 expression, ",i ) ) +
  scale_color_manual( values=c("black","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  xlim(0,16) +
  ylim(0,16) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()

#plot volcano
# label_idx <- ifelse( rownames(qlf$table) %in% c("Sox8","Spib","Tnfrsf11b","Clu","Msln","H2-M2","Nostrin","Pglyrp1",
#                                                 "Ccl9","Ccl6","Ccl20","Tnfaip2","Cr2","Serpin6a","Serpinb1a","Gp2","Cfp","Lyz1","Ctsd"),
#                      rownames(qlf$table), "" )
# color_idx <- ifelse( rownames(qlf$table) %in% c("Sox8","Spib","Tnfrsf11b","Clu","Msln","H2-M2","Nostrin","Pglyrp1",
#                                                 "Ccl9","Ccl6","Ccl20","Tnfaip2","Cr2","Serpin6a","Serpinb1a","Gp2","Cfp","Lyz1","Ctsd"),
#                      "1", "0" )

label_idx <- ifelse( rownames(qlf$table) %in% rownames(topTags(qlf, n=100)), rownames(qlf$table), "" )
color_idx <- ifelse( rownames(qlf$table) %in% rownames(topTags(qlf, n=100)), "1", "0" )

pdf(paste0("figures/",i,"_vs_", j,"_volcano_", date ,".pdf"), height=6, width=9)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("black","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()

fc_cutoff = 1
pval_cutoff = 0.01

de_genes <- abs( qlf$table$logFC ) > fc_cutoff & qlf$table$PValue < pval_cutoff
down_sig <- qlf$table$logFC > fc_cutoff & qlf$table$PValue < pval_cutoff
up_sig <- qlf$table$logFC < -fc_cutoff & qlf$table$PValue < pval_cutoff

write_delim( as.data.frame( rownames(qlf$table)[up_sig] ),
             paste0("de_analysis/de_genes/",i,"_vs_",j,"_UP_de_genes_logfc",fc_cutoff,"_qval",pval_cutoff,".txt"),
             col_names=F,
             delim="\t" )

write_delim( as.data.frame( rownames(qlf$table)[down_sig] ),
             paste0("de_analysis/de_genes/",i,"_vs_",j,"_DOWN_de_genes_logfc",fc_cutoff,"_qval",pval_cutoff,".txt"),
             col_names=F,
             delim="\t" )

#plots with de_gene overlay
# label_idx <- ifelse( rownames(qlf$table) %in% rownames(topTags(qlf, n=50)), "" )

color_idx <- vector( length=nrow(norm_table) )
for (k in 1:nrow(norm_table)) {
  if (rownames(norm_table)[k] %in% rownames(qlf$table)[up_sig]){
    color_idx[k] <- "up"
  } else if (rownames(norm_table)[k] %in% rownames(qlf$table)[down_sig]) {
    color_idx[k] <- "down"
  } else {
    color_idx[k] <- "neutral"
  }
}

pdf(paste0("figures/",i,"_vs_", j,"_exp-exp_up_down_sig_logfc",fc_cutoff,"_qval",pval_cutoff,"_", date ,".pdf"), height=5, width=5)
ggplot( norm_table, aes( x=rowMeans(log2(norm_table[,ctrl_idx])), y=rowMeans(log2(norm_table[,idx])) ) ) +
  geom_point( aes(color=color_idx) ) +
  xlab( paste0( "log2 expression, ",j ) ) +
  ylab( paste0( "log2 expression, ",i ) ) +
  scale_color_manual( values=c("blue","black","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  xlim(0,16) +
  ylim(0,16) +
  ggtitle( paste0(i, " vs ", j) ) +
  annotate( geom="text", label=sum(up_sig), x=2, y=15, color="red" ) +
  annotate( geom="text", label=sum(down_sig), x=15, y=2, color="blue" )
dev.off()

color_idx <- vector( length=nrow(qlf$table) )
for (k in 1:nrow(qlf$table)) {
  if (rownames(qlf$table)[k] %in% rownames(qlf$table)[up_sig]){
    color_idx[k] <- "up"
  } else if (rownames(qlf$table)[k] %in% rownames(qlf$table)[down_sig]) {
    color_idx[k] <- "down"
  } else {
    color_idx[k] <- "neutral"
  }
}

pdf(paste0("figures/",i,"_vs_", j,"_volcano_up_down_sig_logfc",fc_cutoff,"_qval",pval_cutoff,"_", date ,".pdf"), height=4, width=6)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("blue","black","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) ) +
  annotate( geom="text", label=sum(up_sig), x=6, y=5, color="red" ) +
  annotate( geom="text", label=sum(down_sig), x=-6, y=5, color="blue" )
dev.off()

#sig analysis

signature <- read.delim("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures/20210924_Mcell_sig_logfc1_pct10_qval0.01.txt", header=F, sep="\t")
signature <- as.vector(signature$V1)

signature <- read.delim("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de/signatures/20210924_Ionocyte_sig_logfc1_pct10_qval0.01.txt", header=F, sep="\t")
signature <- as.vector(signature$V1)

# signature <- signature[grepl("Rps",signature)==F & grepl("Rpl",signature)==F & grepl("Atp5",signature)==F]

qlf <- qlf.aireko
i = "Mcell.WT.Thy"
j = "Mcell.AireKO.Thy"
idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == i
ctrl_idx <- substr( colnames(count_table), 1, nchar(colnames(count_table)) - 2 ) == j

label_idx <- ifelse( rownames(qlf$table) %in% signature & qlf$table$logFC < -1, rownames(qlf$table), "" )
color_idx <- vector( length=nrow(qlf$table) )
for (k in 1:nrow(qlf$table)) {
  if (rownames(qlf$table)[k] %in% signature){
    color_idx[k] <- "signature+"
  } else {
    color_idx[k] <- "other"
  }
}

pdf(paste0("figures/",i,"_vs_", j,"_volcano_mcell_bespoke_sig_new_lab_", date ,".pdf"), height=4, width=6)
ggplot( qlf$table, aes( x=qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_point( data=qlf$table[signature,], aes( x=qlf$table[signature,"logFC"], 
                                               y=-log10(qlf$table[signature,"PValue"]) ),
              color="purple" ) +
  # geom_text_repel( label=label_idx, color="black", segment.alpha=0.75, segment.color="gray", segment.size=0.25, max.overlaps=25 ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("gray","purple") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  # ggtitle( paste0(i, " vs ", j) ) +
  scale_x_continuous( limits=c(-10,10), oob=squish ) +
  # scale_y_continuous( limits=c(0,5), oob=squish ) +
  annotate( geom="text", label=sum(color_idx=="signature+"&qlf$table$logFC>0), x=6, y=0.5, color="purple" ) +
  annotate( geom="text", label=sum(color_idx=="signature+"&qlf$table$logFC<0), x=-6, y=0.5, color="purple" )
dev.off()


