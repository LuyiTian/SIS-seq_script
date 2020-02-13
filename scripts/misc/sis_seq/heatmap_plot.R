annotation_col

expr = log2(cd[rownames(cd) %in% res_df$gene_names,]+1)

expr = expr[rowSums(expr>0)>10,]

library(Rtsne)
Rtsne_out = Rtsne(expr,dim=1)
Rtsne_out_col = Rtsne(t(expr),dim=1)

rownames(fate_anno_pDC) = fate_anno_pDC$cell_name

pdf("nice_heatmap.pdf")
pheatmap(expr[order(Rtsne_out$Y),],cluster_rows = F,annotation_col = fate_anno_pDC[,c("bias_cDC1", "bias_cDC2", "bias_pDC")],
         show_rownames=F,
         show_colnames=F)
dev.off()

pdf("nice_heatmap_big.pdf",width = 15,height = 25)
pheatmap(expr[order(Rtsne_out$Y),],cluster_rows = F,annotation_col = fate_anno_pDC[,c("bias_cDC1", "bias_cDC2", "bias_pDC")],
         show_rownames=T,
         show_colnames=F)
dev.off()
