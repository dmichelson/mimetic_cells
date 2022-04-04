
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

x.sp <- readRDS( "Rdata/subclustered_postaire_out_2021-02-16.rds" )

args = c( getwd(), "WT1_WT2_KO1_KO2_merged", "2021-06-23", "8" )

library("wesanderson")
library("ggsci")
library("BuenColors")

custom_colormaps = jdb_color_maps
names(custom_colormaps) = seq(1:17)

density_val <- get_density( x.sp@umap[,1], x.sp@umap[,2] )

center = matrix( 0, ncol=2, nrow=13 )
for (i in 1:13) {
    x = mean( x.sp@umap[ x.sp@cluster==i, 1 ] )
    y = mean( x.sp@umap[ x.sp@cluster==i, 2 ] )
    center[i,1] = x
    center[i,2] = y
}

center = as.data.frame( center )

p <- ggplot2::ggplot( as.data.frame(x.sp@umap), ggplot2::aes( x=x.sp@umap[,1], y=x.sp@umap[,2] ) ) +
        ggplot2::geom_point( aes( color = x.sp@cluster, alpha=1 ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        # ggplot2::scale_color_brewer( palette="Paired" ) +
        # ggplot2::scale_color_manual( values = wes_palette( "FantasticFox1", n=12, type="continuous" )  ) +
        ggplot2::scale_color_manual( values = custom_colormaps ) +
        ggplot2::labs(col="gray") +
        ggplot2::ggtitle("Merged MEC Cluster") +
        ggplot2::theme(legend.position="none") +
        ggplot2::geom_point( data = center, mapping = ggplot2::aes_string( x = center[,1], y = center[,2] ), color = "white", size = 6, alpha=0.5, shape="circle" ) +
        ggplot2::geom_text( data = center, mapping = ggplot2::aes_string( x = center[,1], y = center[,2], label = 1:13 ), color = "black", size = 4 )
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_ClusterPlot.pdf" ), height=6, width=6)
print(p)
dev.off()

p <- ggplot2::ggplot( as.data.frame(x.sp@umap), ggplot2::aes( x=x.sp@umap[,1], y=x.sp@umap[,2] ) ) +
        ggplot2::geom_point( aes( color = density_val ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        ggplot2::scale_color_gradientn( colors = jdb_palette("solar_extra")  ) +
        ggplot2::labs(col="density") +
        ggplot2::ggtitle("Merged MEC Cluster") +
        ggplot2::theme()
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_DensityPlot.pdf" ), height=6, width=6)
print(p)
dev.off()

p <- ggplot2::ggplot( as.data.frame(x.sp@umap[x.sp@sample=="WT1" | x.sp@sample=="WT2",]), ggplot2::aes( x=x.sp@umap[x.sp@sample=="WT1"|x.sp@sample=="WT2",1], y=x.sp@umap[x.sp@sample=="WT1"|x.sp@sample=="WT2",2] ) ) +
        ggplot2::geom_point( aes( color = x.sp@cluster[x.sp@sample=="WT1"|x.sp@sample=="WT2"], alpha=1 ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        # ggplot2::scale_color_brewer( palette="Paired" ) +
        # ggplot2::scale_color_manual( values = wes_palette( "FantasticFox1", n=12, type="continuous" )  ) +
        ggplot2::scale_color_manual( values = custom_colormaps ) +
        ggplot2::labs(col="gray") +
        ggplot2::ggtitle("WT MEC Cluster") +
        ggplot2::theme(legend.position="none") +
        ggplot2::geom_point( data = center, mapping = ggplot2::aes_string( x = center[,1], y = center[,2] ), color = "white", size = 6, alpha=0.5, shape="circle" ) +
        ggplot2::geom_text( data = center, mapping = ggplot2::aes_string( x = center[,1], y = center[,2], label = 1:13 ), color = "black", size = 4 )
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_ClusterPlot_WTonly.pdf" ), height=6, width=6)
print(p)
dev.off()

p <- ggplot2::ggplot( as.data.frame(x.sp@umap[x.sp@sample=="KO1" | x.sp@sample=="KO2",]), ggplot2::aes( x=x.sp@umap[x.sp@sample=="KO1"|x.sp@sample=="KO2",1], y=x.sp@umap[x.sp@sample=="KO1"|x.sp@sample=="KO2",2] ) ) +
        ggplot2::geom_point( aes( color = x.sp@cluster[x.sp@sample=="KO1"|x.sp@sample=="KO2"], alpha=1 ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        # ggplot2::scale_color_brewer( palette="Paired" ) +
        # ggplot2::scale_color_manual( values = wes_palette( "FantasticFox1", n=12, type="continuous" )  ) +
        ggplot2::scale_color_manual( values = custom_colormaps ) +
        ggplot2::labs(col="gray") +
        ggplot2::ggtitle("KO MEC Cluster") +
        ggplot2::theme(legend.position="none") +
        ggplot2::geom_point( data = center, mapping = ggplot2::aes_string( x = center[,1], y = center[,2] ), color = "white", size = 6, alpha=0.5, shape="circle" ) +
        ggplot2::geom_text( data = center, mapping = ggplot2::aes_string( x = center[,1], y = center[,2], label = 1:13 ), color = "black", size = 4 )
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_ClusterPlot_KOonly.pdf" ), height=6, width=6)
print(p)
dev.off()

keep <- x.sp@sample=="WT1" | x.sp@sample=="WT2"
density_val <- get_density( x.sp@umap[keep,1], x.sp@umap[keep,2] )

p <- ggplot2::ggplot( as.data.frame(x.sp@umap[x.sp@sample=="WT1" | x.sp@sample=="WT2",]), ggplot2::aes( x=x.sp@umap[x.sp@sample=="WT1"|x.sp@sample=="WT2",1], y=x.sp@umap[x.sp@sample=="WT1"|x.sp@sample=="WT2",2] ) ) +
        ggplot2::geom_point( aes( color = density_val ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        ggplot2::scale_color_gradientn( colors = jdb_palette("solar_extra")  ) +
        ggplot2::labs(col="density") +
        ggplot2::ggtitle("WT MEC Density")
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_DensityPlot_WTonly.pdf" ), height=6, width=6)
print(p)
dev.off()

keep <- x.sp@sample=="KO1" | x.sp@sample=="KO2"
density_val <- get_density( x.sp@umap[keep,1], x.sp@umap[keep,2] )

p <- ggplot2::ggplot( as.data.frame(x.sp@umap[x.sp@sample=="KO1" | x.sp@sample=="KO2",]), ggplot2::aes( x=x.sp@umap[x.sp@sample=="KO1"|x.sp@sample=="KO2",1], y=x.sp@umap[x.sp@sample=="KO1"|x.sp@sample=="KO2",2] ) ) +
        ggplot2::geom_point( aes( color = density_val ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        ggplot2::scale_color_gradientn( colors = jdb_palette("solar_extra")  ) +
        ggplot2::labs(col="density") +
        ggplot2::ggtitle("KO MEC Density")
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_DensityPlot_KOonly.pdf" ), height=6, width=6)
print(p)
dev.off()

#####
#differential density

#####
#differential density plot

# h=c(10,10)

densityWT <- kde2d( x=x.sp@umap[x.sp@sample%in%c("WT1","WT2"),1], 
                    y=x.sp@umap[x.sp@sample%in%c("WT1","WT2"),2],
                    # h=h,
                    n=1000 )
densityKO <- kde2d( x=x.sp@umap[x.sp@sample%in%c("KO1","KO2"),1], 
                    y=x.sp@umap[x.sp@sample%in%c("KO1","KO2"),2],
                    #    h=h,
                    n=1000 )

densityAll <- kde2d( x=x.sp@umap[,1], 
                     y=x.sp@umap[,2],
                    #  h=h,
                    n=1000 )

densityDiff <- densityWT$z - densityKO$z

xy = xy.coords( x=x.sp@umap[,1],
                y=x.sp@umap[,2] )
select = is.finite(xy$x) & is.finite(xy$y)
x = cbind(xy$x, xy$y)[select, ]
mkBreaks = function(u) u - diff(range(u))/(length(u) - 1)/2
xbin = cut(x[, 1], mkBreaks(densityAll$x), labels = FALSE)
ybin = cut(x[, 2], mkBreaks(densityAll$y), labels = FALSE)
dens = densityDiff[cbind(xbin, ybin)]
dens[is.na(dens)] = 0

pdf( "figures/2021-02-16_wt_ko_scatac_diff_density.pdf", height=6, width=7 )
ggplot( as.data.frame( x.sp@umap ),
        aes( x=x.sp@umap[,1],
             y=x.sp@umap[,2] ) ) +
  geom_point( aes( color=squish(dens, range=c(-0.03,0.03)) ), size=0.5 ) +
  theme_bw() +
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
  scale_color_gradient2( low="purple4", mid="navajowhite2", high="green4" ) +
  labs( color="Density differential\n(WT - KO)" ) +
  xlab("UMAP-1") +
  ylab("UMAP-2")
dev.off()

#####
#only 1 and 7 
dat <- x.sp[x.sp@cluster%in%c(1,7),,]
p <- ggplot2::ggplot( as.data.frame(dat@umap), ggplot2::aes( x=dat@umap[,1], y=dat@umap[,2] ) ) +
        ggplot2::geom_point( aes( color = dat@cluster, alpha=1 ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        # ggplot2::scale_color_brewer( palette="Paired" ) +
        # ggplot2::scale_color_manual( values = wes_palette( "FantasticFox1", n=12, type="continuous" )  ) +
        ggplot2::scale_color_manual( values = c("#00441B","#FFA300") ) +
        ggplot2::labs(col="gray") +
        ggplot2::ggtitle("Merged MEC Cluster") +
        ggplot2::theme(legend.position="none") +
        ggplot2::geom_point( data = center[c(1,7),], mapping = ggplot2::aes_string( x = center[c(1,7),1], y = center[c(1,7),2] ), color = "white", size = 6, alpha=0.5, shape="circle" ) +
        ggplot2::geom_text( data = center[c(1,7),], mapping = ggplot2::aes_string( x = center[c(1,7),1], y = center[c(1,7),2], label = c(1,7) ), color = "black", size = 4 )
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_ClusterPlot_1v7only_noBackground.pdf" ), height=6, width=6)
print(p)
dev.off()

dat <- x.sp
dat@cluster[dat@cluster%in%c("1","7")==F] <- "8"
p <- ggplot2::ggplot( as.data.frame(dat@umap), ggplot2::aes( x=dat@umap[,1], y=dat@umap[,2] ) ) +
        ggplot2::geom_point( aes( color = dat@cluster, alpha=1 ) ) +
        ggplot2::geom_point( data=as.data.frame(dat@umap[dat@cluster%in%c("1","7"),]), aes( x=dat@umap[dat@cluster%in%c("1","7"),1], y=dat@umap[dat@cluster%in%c("1","7"),2],color = dat@cluster[dat@cluster%in%c("1","7")], alpha=1 ) ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        # ggplot2::scale_color_brewer( palette="Paired" ) +
        # ggplot2::scale_color_manual( values = wes_palette( "FantasticFox1", n=12, type="continuous" )  ) +
        ggplot2::scale_color_manual( values = c("#00441B","#FFA300","gray") ) +
        ggplot2::labs(col="gray") +
        ggplot2::ggtitle("Merged MEC Cluster") +
        ggplot2::theme(legend.position="none") +
        ggplot2::geom_point( data = center[c(1,7),], mapping = ggplot2::aes_string( x = center[c(1,7),1], y = center[c(1,7),2] ), color = "white", size = 6, alpha=0.5, shape="circle" ) +
        ggplot2::geom_text( data = center[c(1,7),], mapping = ggplot2::aes_string( x = center[c(1,7),1], y = center[c(1,7),2], label = c(1,7) ), color = "black", size = 4 )
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_ClusterPlot_1v7only_withBackground.pdf" ), height=6, width=6)
print(p)
dev.off()

dat <- x.sp
dat@cluster[dat@cluster%in%c("8")==F] <- "1"
p <- ggplot2::ggplot( as.data.frame(dat@umap), ggplot2::aes( x=dat@umap[,1], y=dat@umap[,2] ) ) +
        ggplot2::geom_point( aes( color = dat@cluster, alpha=1 ) ) +
        ggplot2::geom_point( data=as.data.frame(dat@umap[dat@cluster==8,]), aes( x=dat@umap[dat@cluster==8,1], y=dat@umap[dat@cluster==8,2] ), color="#C390D4" ) +
        ggplot2::theme_classic() +
        ggplot2::xlab("UMAP-1") +
        ggplot2::ylab("UMAP-2") +
        # ggplot2::scale_color_brewer( palette="Paired" ) +
        # ggplot2::scale_color_manual( values = wes_palette( "FantasticFox1", n=12, type="continuous" )  ) +
        ggplot2::scale_color_manual( values = c("gray","#C390D4") ) +
        ggplot2::labs(col="gray") +
        ggplot2::ggtitle("Merged MEC Cluster") +
        ggplot2::theme(legend.position="none") +
        ggplot2::geom_point( data = center[c(8),], mapping = ggplot2::aes_string( x = center[c(8),1], y = center[c(8),2] ), color = "white", size = 6, alpha=0.5, shape="circle" ) +
        ggplot2::geom_text( data = center[c(8),], mapping = ggplot2::aes_string( x = center[c(8),1], y = center[c(8),2], label = c(8) ), color = "black", size = 4 )
pdf(paste0( "figures/", args[2], "_", args[3],"_umap_ClusterPlot_8only_withBackground.pdf" ), height=6, width=6)
print(p)
dev.off()

#####
#hclust

# hierarchical clustering
# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
	SnapATAC::colMeans(x.sp[x,], mat="bmat")
} )

# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls))/2)), method="mcquitty");

pdf(paste( "figures/", args[2],"_",args[3], "_hclust_viz.pdf", sep="" ), height=5, width=6)
# plotViz(
#     obj=x.sp,
#     method="umap", 
#     main="Merged MEC Cluster",
#     point.color=x.sp@cluster, 
#     point.size=1, 
#     point.shape=19, 
#     point.alpha=0.8, 
#     text.add=TRUE,
#     text.size=1.5,
#     text.color="black",
#     text.halo.add=TRUE,
#     text.halo.color="white",
#     text.halo.width=0.2,
#     down.sample=10000,
#     legend.add=FALSE
# 	)
plot(hc, hang=-1, xlab="")
plot(as.phylo(hc), cex = 0.6, label.offset = 0.5)
plot(as.phylo(hc), type = "cladogram", cex = 0.6, label.offset = 0.5)
plot(as.phylo(hc), type = "fan")
plot(as.phylo(hc),type="unrooted", cex=0.6, no.margin=T)
dev.off()

dist_mat = cor(t(do.call(rbind, ensemble.ls)))

breaksList = seq(0.2,1, by=0.05)
pdf(paste( "figures/", args[2],"_",args[3], "_distmat_viz.pdf", sep="" ), height=4.75, width=4)
pheatmap(dist_mat,
        cluster_rows=hc,
        cluster_cols=hc,
        treeheight_row=0,
        treeheight_col=100,
        breaks=breaksList,
        color=colorRampPalette(brewer.pal(name="Reds",n=9))(length(breaksList)),
        border_color=NA
)
dev.off()