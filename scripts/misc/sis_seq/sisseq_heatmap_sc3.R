
library(SC3)
library(pheatmap)

sce_sc3 = sc3(sce_sc3[rownames(sce_sc3) %in% sel_genes,],
              ks=3:8,pct_dropout_min = 15, pct_dropout_max = 95,k_estimator=T,rand_seed=66666,biology=T)

sce_sc3_hi_var = sc3(sce_sc3_hi_var,
              ks=3:8,pct_dropout_min = 15, pct_dropout_max = 95,k_estimator=T,rand_seed=66666,biology=T)

organise_marker_genes.remake <- function(object, k, p_val, auroc) {
  dat <- rowData(object)[, c(paste0("sc3_", k, "_markers_clusts"), paste0("sc3_", k, 
                                                                          "_markers_auroc"), paste0("sc3_", k, "_markers_padj"), "feature_symbol")]
  dat <- dat[dat[, paste0("sc3_", k, "_markers_padj")] < p_val & !is.na(dat[, paste0("sc3_", 
                                                                                     k, "_markers_padj")]), ]
  dat <- dat[dat[, paste0("sc3_", k, "_markers_auroc")] > auroc, ]
  
  d <- NULL
  
  for (i in sort(unique(dat[, paste0("sc3_", k, "_markers_clusts")]))) {
    tmp <- dat[dat[, paste0("sc3_", k, "_markers_clusts")] == i, ]
    tmp <- tmp[order(tmp[, paste0("sc3_", k, "_markers_auroc")], decreasing = TRUE), ]
    d <- rbind(d, tmp)
  }
  
  return(d)
}


markers_for_heatmap.remake <- function(markers) {
  res <- NULL
  for (i in unique(markers[, 1])) {
    tmp <- markers[markers[, 1] == i, ]
    if (nrow(tmp) > 10) {
      res <- rbind(res, tmp[1:10, ])
    } else {
      res <- rbind(res, tmp)
    }
  }
  
  return(res)
}

colnames(sce_sc3) = sce_sc3$cell_name

sc3_plot_markers(sce_sc3,k=8,show_pdata=c("bias_cDC1","bias_cDC2","bias_pDC"),p.val=0.05,auroc=0.6)

sc3_plot_markers(sce_sc3_hi_var,k=8,show_pdata=c("bias_cDC1","bias_cDC2","bias_pDC"),p.val=0.05,auroc=0.6)

markers = organise_marker_genes.remake(sce_sc3,k=4,p_val=1,auroc=0.1)

expr = logcounts(sce_sc3)[markers$feature_symbol,]

library(Rtsne)
set.seed(66666)
tsne_out_gene = Rtsne((expr),dim=1)
tsne_out_cell = Rtsne(t(expr),dim=1)

expr_tsne = expr[order(tsne_out_gene$Y),order(tsne_out_cell$Y)]

col_anno = as.data.frame(colData(sce_sc3)[,c("bias_cDC1","bias_cDC2","bias_pDC")])
#rownames(col_anno) = rownames(colData(sce_sc3))

ann_colors = list(
  bias_cDC1 = c("#FFFFFF", "#FF7272"),
  bias_cDC2 = c("#FFFFFF", "#74FF72"),
  bias_pDC = c("#FFFFFF", "#7276FF")
)
k=7
hc <- metadata(sce_sc3)$sc3$consensus[[as.character(k)]]$hc

pdf("sc3_sis_seq_heatmap.pdf",width = 8,height = 16)
pheatmap(expr, annotation_col=col_anno,annotation_colors=ann_colors, show_colnames = FALSE)
dev.off()


pheatmap(expr, annotation_col=col_anno,annotation_colors=ann_colors, cluster_cols = hc,show_colnames = FALSE)

pheatmap(expr_tsne, cluster_rows=F,cluster_cols=F,annotation_col=col_anno,annotation_colors=ann_colors ,show_colnames = FALSE)




uni_genes = intersect(rownames(sce_TabulaMuris_all),rownames(expr))


pdf("sc3_sis_seq_heatmap_TM_intersect.pdf",width = 8,height = 16)
pht = pheatmap(expr[uni_genes,], annotation_col=col_anno,annotation_colors=ann_colors, show_colnames = FALSE)
dev.off()

colnames(sce_TabulaMuris_all) = sce_TabulaMuris_all$cell
tmp = as.data.frame(colData(sce_TabulaMuris_all)[,c("cell_ontology_class","total_features")])
rownames(tmp) = rownames(colData(sce_TabulaMuris_all))


pdf("sc3_sis_seq_TabulaMuris_heatmap.pdf",width = 20,height = 16)
pheatmap(logcounts(sce_TabulaMuris_all)[uni_genes, ],cluster_rows=pht$tree_row,annotation_col=tmp)
dev.off()
