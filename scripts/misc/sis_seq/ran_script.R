selected_genes_3fate <- read_csv("Dropbox/research/sis_seq/selected_genes_3fate.csv")

head(selected_genes_3fate)
dim(selected_genes_3fate)
any(duplicated(selected_genes_3fate$gene_names))
dim(selected_genes_3fate[!duplicated(selected_genes_3fate$gene_names),])
library(readr)
selected_gRNA_broadgpp_brie_library_457X4genes <- read_csv("Dropbox/research/sis_seq/selected_gRNA_broadgpp.brie.library_457X4genes.csv")
selected_gRNA_broadgpp_brie_library_457X4genes <- read.csv("~/Downloads/selected_gRNA_broadgpp.brie.library_457X4genes.csv", stringsAsFactors=FALSE)
selected_gRNA_broadgpp_brie_library_457X4genes$Target.Gene.Symbol[!(selected_gRNA_broadgpp_brie_library_457X4genes$Target.Gene.Symbol %in% selected_genes_3fate$gene_names)]
unique(selected_gRNA_broadgpp_brie_library_457X4genes$Target.Gene.Symbol[!(selected_gRNA_broadgpp_brie_library_457X4genes$Target.Gene.Symbol %in% selected_genes_3fate$gene_names)])

cDC1_genes <- read_csv("Dropbox/research/sis_seq/cDC1_genes.csv",
col_names = FALSE)


pDC_genes <- read_csv("Dropbox/research/sis_seq/pDC_genes.csv",
col_names = FALSE)


cDC2_genes <- read_csv("Dropbox/research/sis_seq/cDC2_genes.csv",col_names = FALSE)


sis_seq_all = unique(c(cDC1_genes$X1,cDC2_genes$X1,pDC_genes$X1))

unique(selected_gRNA_broadgpp_brie_library_457X4genes$Target.Gene.Symbol[!(selected_gRNA_broadgpp_brie_library_457X4genes$Target.Gene.Symbol %in% sis_seq_all)])

known_list = c("Cbfa2t3","Mycl","Id2","Batf3","Irf8","Runx1","Zbtb46","Runx2","Tcf4","Zeb2","Gfi1","Ltbr","Notch2","Klf4","Spib")

clu1 = srt_sub.markers[srt_sub.markers$cluster==1 & srt_sub.markers$avg_diff>0,] %>% top_n(n=60,wt=avg_diff)
clu3 = srt_sub.markers[srt_sub.markers$cluster==3 & srt_sub.markers$avg_diff>0,] %>% top_n(n=60,wt=avg_diff)
clu5 = srt_sub.markers[srt_sub.markers$cluster==5 & srt_sub.markers$avg_diff>0,] %>% top_n(n=60,wt=avg_diff)

write.csv(Reduce(intersect,list(clu1$gene,clu3$gene,clu5$gene)),quote=FALSE,row.names = FALSE)

write.csv(clu1$gene,quote=FALSE,row.names = FALSE)

write.csv(clu3$gene,quote=FALSE,row.names = FALSE)


write.csv(clu5$gene,quote=FALSE,row.names = FALSE)



tmp_wrt = srt_sub.markers[srt_sub.markers$p_val_adj<0.01 & srt_sub.markers$avg_diff>0,] %>% group_by(cluster) %>% top_n(n=60,wt=avg_diff)
write.csv(tmp_wrt,file="supp_table_top60_marker_per_clusters.csv",row.names = FALSE)






