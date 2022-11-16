library(ggplot2)
library(dplyr)

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/text_outputs/2021-09-16_sig_de")
date="20210924"

sig_list <- list()

for (i in list.files()) {
  if (substr(i,nchar(i)-3,nchar(i)) ==".txt") {
    temp_table <- read.delim(i,header=T,sep="\t")
    temp_idx <- temp_table$avg_log2FC > 1 & temp_table$pct.1 > 0.1 & temp_table$p_val_adj < 0.01
    sig_list[[substr(i,25,nchar(i)-15)]] <- temp_table[temp_idx,1]
  }
}

sig_list_unique <- list()

for (i in names(sig_list)) {
  temp_sig <- sig_list[[i]]
  for (j in names(sig_list)) {
    if (i==j) { next }
    else {
      temp_idx <- temp_sig %in% sig_list[[j]]
      temp_sig <- temp_sig[!temp_idx]
    }
    sig_list_unique[[i]] <- temp_sig
  }
}

for (i in names(sig_list_unique)) {
  write_delim( as.data.frame(sig_list_unique[[i]]), paste0("signatures/",date,"_",i,"_sig_logfc1_pct10_qval0.01.txt"), col_names=F, delim="\t" )
}
