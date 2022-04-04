calculate_cluster_fraction <- function( n_hash, names_hash, n_cluster, names_cluster, n_geno, names_geno, 
                                        hash_path, cluster_path ) {
  
  cluster_frac <- matrix( nrow=n_cluster, ncol=n_hash )
  colnames(cluster_frac) <- names_hash
  rownames(cluster_frac) <- names_cluster
  
  for (i in names_cluster) {
    for (j in names_hash) {
      cluster_frac[i,j] <- sum( hash_path==j & cluster_path==i ) / sum( hash_path==j )
    }
  }
  
  temp <- vector()
  for (i in ncol(cluster_frac)) {
    temp <- c(temp, cluster_frac[,i])
  }
  
  cluster_frac <- data.frame( frac=as.numeric(cluster_frac), 
                              geno=factor(rep(substr(names_hash, 1, nchar(as.vector(names_hash))-3), each=n_cluster), 
                                          levels=names_geno, ordered=T),
                              hash=factor(rep(names_hash, each=n_cluster), levels=names_hash, ordered=T), 
                              cluster=factor(names_cluster, ordered=T) )

  mean_cluster_frac <- matrix(ncol=3, nrow=n_cluster*n_geno)
  count=0
  colnames(mean_cluster_frac) <- c("frac","genotype","cluster")
  for (i in names_geno) {
    for (j in names_cluster ){
      count=count+1
      mean_cluster_frac[count,"genotype"] = i
      mean_cluster_frac[count,"cluster"] = j
      mean_cluster_frac[count,"frac"] = mean( cluster_frac[grepl(i,cluster_frac$hash) & cluster_frac$cluster==j,"frac"] )
    }
  }
  
  mean_cluster_frac <- data.frame( "frac"=as.numeric(mean_cluster_frac[,1]),
                                   "genotype"=factor(mean_cluster_frac[,2], levels=names_geno, ordered=T),
                                   "cluster"=factor(mean_cluster_frac[,3],levels=names_cluster, ordered=T) )
  
  return(list(cluster_frac, mean_cluster_frac))
  
}
