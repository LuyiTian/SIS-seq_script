





var.fit <- trendVar(sce, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce, var.fit)
lvg.out <- var.out[order(var.out$bio)[1:3000], ]
lvg.out = lvg.out[lvg.out$mean>3.5,]

pdf("select_lvg.pdf")
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
points(var.out$mean[rownames(var.out) %in% rownames(lvg.out)], var.out$total[rownames(var.out) %in% rownames(lvg.out)], col="red", pch=16)
dev.off()

comm_gene = intersect(rownames(lvg.out), rownames(gene_expr_t))
SISseq_DC_expr = logcounts(sce)[comm_gene, ]
HSC_DC_expr = gene_expr_t[comm_gene, ]

corr_mat = apply(HSC_DC_expr,2,function(x){cor(x,SISseq_DC_expr)})


colSums(HSC_DC_expr)[order(colSums(HSC_DC_expr),decreasing = T)][1:50]


pdf("big_HSC_DC_heatmap.pdf",height = 15,width = 15)
pheatmap(HSC_DC_expr[,grepl("H",colnames(HSC_DC_expr))],annotation_col=col_anno)
dev.off()



cell_cycle_genes = get_genes_by_GO(returns = "external_gene_name",dataset = "mmusculus_gene_ensembl", go = "GO:0007049")

HSC_cc = rownames(HSC_DC_expr)[rownames(HSC_DC_expr) %in% cell_cycle_genes]

HSC_gene_corr = cor(t(HSC_DC_expr))

HSC_cc_ext = c()

for(g in HSC_cc){
  tmp = HSC_gene_corr[g,]
  tmp = tmp[order(abs(tmp),decreasing = T)[1:4]]
  HSC_cc_ext = c(HSC_cc_ext,names(tmp))
}

HSC_cc_ext = unique(HSC_cc_ext)
HSC_cc = unique(c(HSC_cc, HSC_cc_ext))

col_anno = data.frame(cell_type=coordinates_gene_counts_flow_cytometry$group)
rownames(col_anno) = colnames(HSC_DC_expr)
pdf("big_HSC_prog_DC_heatmap_cc_rm.pdf",height = 15,width = 15)
pheatmap(HSC_DC_expr[!(rownames(HSC_DC_expr) %in% HSC_cc),],annotation_col=col_anno)
dev.off()

tmp_expr = HSC_DC_expr[!(rownames(HSC_DC_expr) %in% HSC_cc),grepl("H",colnames(HSC_DC_expr))]
t_out = Rtsne(t(tmp_expr),dim=1)

pdf("big_HSC_DC_heatmap_cc_rm.pdf",height = 15,width = 15)
pheatmap(tmp_expr[,order(t_out$Y)],annotation_col=col_anno,show_colnames=F,cluster_cols=F)
dev.off()


