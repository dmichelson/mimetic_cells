library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(pheatmap)

setwd("C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_RNA-seq/2020-08-21_mcell_meclo_rnaseq/Batch096_Michelson_Mouse")

date = "2021-04-05"

norm_table <- read.delim( "gene_expression_1_adjust.gct", sep="\t", header=T, skip=2, row.names=1 )
norm_table <- norm_table[,c(5:7,14:19)]
count_table <- read.delim( "genes_count_table_QC_passed.csv", sep=",", header=T, row.names=1 )
count_table <- count_table[,c(4:6,13:18)]

sample_names <- colnames(norm_table)
 
setwd("meclo")

#####
#look at counts data

#look at counts per sample
pdf( paste0("figures/counts_per_sample_",date,".pdf"), height=3, width=5 )
hist(log10(colSums(count_table)), breaks=3, xlim=range(0,10))
abline( v=3.3010, col="red" )
dev.off()

#####
#look at normalized data

pdf( paste0("figures/genes_per_sample_",date,".pdf"), height=3, width=5 )
hist(colSums( norm_table > 1 ), breaks=3, xlim=range(0,50000))
dev.off()

# remove genes with fewer than two samples with at least expression 10
pdf( paste0("figures/detection_per_gene_",date,".pdf"), height=3, width=5 )
hist(rowSums( norm_table > 10 ), breaks=10)
abline(v=2, col="red")
dev.off()

idx <- rowSums( norm_table > 10 ) < 2
norm_table <- norm_table[!idx,]

#####
# calculate basic stats for samples and genes
means <- rowMeans(norm_table)
vars <- apply(norm_table,1,var)
cvs <- apply(norm_table,1,sd)/rowMeans(norm_table)

pdf( paste0("figures/cv_vs_mean_",date,".pdf"), height=5, width=5 )
plot(log(means),log(cvs))
dev.off()


#####
#look at sample clustering
sample_cor_matrix <- cor( norm_table, method="pearson" )
hclust_samples <- hclust( as.dist( 1-sample_cor_matrix ), method="ward.D2" )

gene_cor_matrix <- cor( t( norm_table ), method = "spearman" )
hclust_genes <- hclust( as.dist( 1-gene_cor_matrix ), method="ward.D2" )

scaled_norm_table <- log2( norm_table )
scaled_norm_table <- scaled_norm_table - rowMeans(scaled_norm_table)

annotation_names = as.data.frame(sample_names)
rownames( annotation_names ) = rownames( sample_cor_matrix )

breaksList <- seq(-1,1,by=0.01)

pdf( paste0("figures/samples_distance_heatmap_unclustered_", date, ".pdf"), height=5, width=7.5)
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
          breaks=seq(0.2,1,by=0.01), 
          color=colorRampPalette((brewer.pal(name="Reds",n=9)))(80),
          border_color = NA
)
dev.off()

pdf( paste0("figures/samples_distance_heatmap_clustered_", date, ".pdf"), height=5, width=7.5)
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
          breaks=seq(0.2,1,by=0.01), 
          color=colorRampPalette((brewer.pal(name="Reds",n=9)))(80),
          border_color = NA
)
dev.off()


#####
#pca analysis
library(pcaMethods)
library(ggrepel)

pca_object <- pca( t(norm_table), method="svd", nPcs=6, scale="uv", center=T )

pdf( paste0("figures/elbow_plot_", date, ".pdf" ), height=3, width=3 )
plot( pca_object@R2, xlab="PC#", ylab="R^2", main="Elbow plot" )
lines( pca_object@R2 )
dev.off()

label_idx <- c("Aire.MEC.1","Aire.MEC.2","Aire.MEC.3",
               "pre-Aire.MEC.1","pre-Aire.MEC.2","pre-Aire.MEC.3",
               "post-Aire.MEC.1","post-Aire.MEC.2","post-Aire.MEC.3")

pdf( paste0("figures/PCA_plot_", date, ".pdf" ), height=3, width=4 )
for (i in 1:5) {
  temp_plot <- ggplot( as.data.frame(pca_object@scores), aes( x=pca_object@scores[,i],y=pca_object@scores[,i+1] ) ) +
    geom_point( aes( fill=factor( substr(sample_names,1,nchar(sample_names)-2 ) ) ), size=3, alpha=1, shape=21 ) +
    geom_text_repel( label=label_idx ) +
    theme_classic() +
    theme( legend.position="none" ) +
    xlab( paste0("PC",i, " (",pca_object@R2[i]*100,"%)") ) +
    ylab( paste0("PC",i+1, " (",pca_object@R2[i+1]*100,"%)") ) +
    theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5) ) +
    scale_fill_manual( values=c("black","blue","orange") )
  print( temp_plot )
}
dev.off()

#find loadings for PCs
pca_genes <- matrix( ncol=10, nrow=50 )
colnames( pca_genes ) = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
for (i in 1:6) {
  temp <- sort( pca_object@loadings[, paste0("PC",i) ], decreasing=T )
  pca_genes[,i] = names( temp[1:50] )
}

write_delim( as.data.frame( pca_genes ), paste0( "de_analysis/top_pca_genes_per_component_", date, ".txt" ), delim="\t", col_names=T )

pca_genes <- matrix( ncol=10, nrow=50 )
colnames( pca_genes ) = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
for (i in 1:6) {
  temp <- sort( pca_object@loadings[, paste0("PC",i) ], decreasing=F )
  pca_genes[,i] = names( temp[1:50] )
}

write_delim( as.data.frame( pca_genes ), paste0( "de_analysis/bottom_pca_genes_per_component_", date, ".txt" ), delim="\t", col_names=T )

#####
# de analysis
i = "PDPNn.CD104n.MEClo.Thy"
j = "PDPNp.CD104p.MEClo.Thy"

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

write_delim( as.data.frame( cbind( rownames(qlf$table),qlf$table ) ), paste0("de_analysis/",i,"_vs_",j,"_de_table.txt"), col_names=T, delim="\t" )

#plot volcano
x <- c("Hnf4a","Hnf4g","Spdef","Foxa1","Foxa2","Foxa3","Grhl1","Grhl3",
       "Spib","Sox8","Foxj1","Myog","Gata3","Pou2f3","Foxi1","Foxi2","Ptf1a")

label_idx <- ifelse( rownames(qlf$table) %in% x, rownames(qlf$table), "" )
color_idx <- ifelse( rownames(qlf$table) %in% x, "1", "0" )

pdf( paste0("figures/",i,"_vs_",j,"_TF_volcano_",date,".pdf"), height=3, width=4)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$q_val) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_point( data=qlf$table[x,], aes( x=-qlf$table[x,"logFC"], y=-log10(qlf$table[x,"q_val"]) ), color="darkred" ) +
  geom_label_repel( label=label_idx, color="red", force=0.8 ) +
  xlab( paste0("log2 fold change, post-Aire vs pre-Aire MEC") ) +
  ylab("-log10 adjusted P-value") +
  scale_color_manual( values=c("black","red") ) +
  theme_classic() +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.5), legend.position="none" ) +
  ggtitle( paste0( "Post-Aire vs pre-Aire MEC") )
dev.off()

fc_cutoff = 1
qval_cutoff = 0.05

de_genes <- abs( qlf$table$logFC ) > fc_cutoff & qlf$table$q_val < qval_cutoff
down_sig <- qlf$table$logFC > fc_cutoff & qlf$table$q_val < qval_cutoff
up_sig <- qlf$table$logFC < -fc_cutoff & qlf$table$q_val < qval_cutoff

write_delim( as.data.frame( rownames(qlf$table)[up_sig] ),
             paste0("de_analysis/de_genes/",i,"_vs_",j,"_UP_de_genes_logfc",fc_cutoff,"_qval",qval_cutoff,".txt"),
             col_names=F,
             delim="\t" )

write_delim( as.data.frame( rownames(qlf$table)[down_sig] ),
             paste0("de_analysis/de_genes/",i,"_vs_",j,"_DOWN_de_genes_logfc",fc_cutoff,"_qval",qval_cutoff,".txt"),
             col_names=F,
             delim="\t" )


#de volcano plot
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

pdf(paste0("figures/",i,"_vs_", j,"_volcano_up_down_sig_logfc",fc_cutoff,"_qval",qval_cutoff,"_", date ,".pdf"), height=4, width=6)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$q_val) ) ) +
  geom_point( aes(color=color_idx) ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 adjusted P-value") +
  scale_color_manual( values=c("blue","gray","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) ) +
  xlim(-10,10) +
  annotate( geom="text", label=sum(up_sig), x=10, y=7, color="red" ) +
  annotate( geom="text", label=sum(down_sig), x=-10, y=7, color="blue" )
dev.off()

#signature analysis
aig_table <- read.delim( "AIG_gene_matched.txt", sep="\t", header=T )
ang_table <- read.delim( "ANG_gene_matched.txt", sep="\t", header=T )

aig_index <- rownames(qlf$table)%in%aig_table[,1]
ang_index <- rownames(qlf$table)%in%ang_table[,1]

aig_repressed <- sum( aig_index==T & qlf$table$logFC > 0 )
aig_induced <- sum( aig_index==T & qlf$table$logFC < 0 )

ang_repressed <- sum( ang_index==T & qlf$table$logFC > 0 )
ang_induced <- sum( ang_index==T & qlf$table$logFC < 0 )

total_repressed <- sum( qlf$table$logFC > 0 )
total_induced <- sum( qlf$table$logFC < 0 )

color_idx <- vector(length=nrow(qlf$table))
for (k in 1:nrow(qlf$table)) {
  if ( aig_index[k]==T ) {
    color_idx[k] <- "aig"
  } else if ( ang_index[k]==T ) {
    color_idx[k] <- "ang"
  } else {
    color_idx[k] <- "other"
  }
}

pdf(paste0("figures/",i,"_vs_", j,"_volcano_AIG_ANG_", date ,".pdf"), height=3, width=4.5)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$q_val) ) ) +
  geom_point( aes(color=color_idx) ) +
  xlab( paste0( "log2 fold change, post-Aire vs pre-Aire MEC" )) +
  ylab("-log10 adjusted P-value") +
  scale_color_manual( values=c("red","blue","black") ) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  annotate( geom="text", label=aig_induced, x=10, y=7.5, color="red" ) +
  annotate( geom="text", label=aig_repressed, x=-8, y=7.5, color="red" ) +
  annotate( geom="text", label=ang_induced, x=10, y=7, color="blue" ) +
  annotate( geom="text", label=ang_repressed, x=-8, y=7, color="blue" ) +
  annotate( geom="text", label=total_induced, x=10, y=6.5, color="black" ) +
  annotate( geom="text", label=total_repressed, x=-8, y=6.5, color="black" ) +
  labs( color="" ) +
  xlim(-9,11) +
  ggtitle( paste0( i, " vs ", j ) )
dev.off()
