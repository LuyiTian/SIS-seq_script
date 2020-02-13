# heatmap

library(pheatmap)
tt_cnt = edgeR::cpm(allcounts1$counts)
gm_cnt = log2(tt_cnt[rownames(allcounts1$counts) %in% unique(c(pDC_genes, cDC1_genes, cDC2_genes)), ]+1)
#gm_cnt = log2(tt_cnt[unique(rbind(pDC_sel,cDC1_sel, cDC2_sel)$gene_names), ]+1)

gm_cnt = gm_cnt[rowMeans(gm_cnt)<6 & rowSums(gm_cnt>0)<70,]


tmp_col_anno = fate_anno_pDC[,c("bias_cDC1", "bias_cDC2", "bias_pDC")]
rownames(tmp_col_anno) = fate_anno_pDC$cell_name

pdf("pheatmap_all.pdf")
pheatmap(gm_cnt, annotation_col=tmp_col_anno)
dev.off()

library(Rmagic)

MAGIC_data = magic(t(gm_cnt), k=4,t=3)
rownames(MAGIC_data) = colnames(gm_cnt)
MAGIC_data = scale(MAGIC_data,center = TRUE, scale = TRUE)
#MAGIC_data = MAGIC_data[,colSums(MAGIC_data) > 0]
MAGIC_data[is.na(MAGIC_data)]=0
pdf("pheatmap_all_MAGIC.pdf")

pheatmap(t(as.matrix(MAGIC_data)),annotation_col=tmp_col_anno)

dev.off()

tsne_magic_out_cell = Rtsne::Rtsne(MAGIC_data,dim=2,perplexity=5)
tsne_magic_out_gene = Rtsne::Rtsne(t(MAGIC_data),dim=1,perplexity=5,check_duplicates=F)
pheatmap(t(as.matrix(MAGIC_data))[order(tsne_magic_out_gene$Y),],annotation_col=tmp_col_anno,cluster_rows = F,cluster_cols = hclust(dist(tsne_magic_out_cell$Y)))

library(scran)

mnn_result = mnnCorrect(gm_cnt[,fate_anno_pDC$cell_name[fate_anno_pDC$batch=="JS86"]],
           gm_cnt[,fate_anno_pDC$cell_name[fate_anno_pDC$batch=="JT27"]],
           gm_cnt[,fate_anno_pDC$cell_name[fate_anno_pDC$batch=="JS91"]],
           gm_cnt[,fate_anno_pDC$cell_name[fate_anno_pDC$batch=="JS96"]],k=5,sigma=0.5)

mnn_cnt = Reduce(cbind,mnn_result$corrected)
colnames(mnn_cnt) = colnames(gm_cnt)
mnn_cnt[mnn_cnt<0] = 0
pdf("pheatmap_all_MNN.pdf")

pheatmap(mnn_cnt, annotation_col=tmp_col_anno)

dev.off()





gm_cnt = log2(tt_cnt[rownames(allcounts1$counts) %in% unique(c(pDC_genes, cDC1_genes, cDC2_genes)), ]+1)



library(SC3)

sce_sc3 = SingleCellExperiment(assay=list(counts=allcounts1$counts))
colData(sce_sc3) = DataFrame(fate_anno_pDC)

sce_sc3 <- computeSumFactors(sce_sc3,sizes=seq(5, 50, 5))
sce_sc3 <- normalize(sce_sc3)

var.fit <- trendVar(sce_sc3, method="loess", use.spikes=FALSE)
var.out <- decomposeVar(sce_sc3, var.fit)
hvg.out <- var.out[order(var.out$bio, decreasing=TRUE)[1:1500], ]
hvg.out = hvg.out[hvg.out$mean>1 & hvg.out$mean<9,]

plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression",
     ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
points(var.out$mean[rownames(var.out) %in% rownames(hvg.out)], var.out$total[rownames(var.out) %in% rownames(hvg.out)], col="red", pch=16)


rowData(sce_sc3)$feature_symbol = rownames(sce_sc3)

sce_sc3 = sc3(sce_sc3[rownames(sce_sc3) %in% rownames(hvg.out),],
    ks=3:8,pct_dropout_min = 10, pct_dropout_max = 95,k_estimator=T,rand_seed=66666,biolog=T)


######
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
######


sc3_plot_markers(sce_sc3,k=7,show_pdata=c("bias_cDC1","bias_cDC2","bias_pDC"),p.val=0.05,auroc=0.6)



#####

library(DrImpute)
imp <- DrImpute(logcounts(sce_sc3))
assay(sce_sc3, "scran_DrImpute") <- imp

pheatmap(assay(sce_sc3, "scran_DrImpute")[rownames(sce_sc3) %in% sel_genes,])

######



######
com_gene = intersect(rownames(gm_cnt),rownames(t(MAGIC_data)))
ct1 = gm_cnt[com_gene,]
ct2 = t(MAGIC_data)[com_gene,]
sce_sc3_magic = SingleCellExperiment(assay=list(counts=ct1,logcounts=ct2))
sce_sc3_magic <- computeSumFactors(sce_sc3_magic,sizes=seq(5, 50, 5))
sce_sc3_magic <- normalize(sce_sc3_magic,return_log=F)

colData(sce_sc3_magic) = DataFrame(fate_anno_pDC)
rowData(sce_sc3_magic)$feature_symbol = rownames(sce_sc3_magic)

sce_sc3_magic = sc3(sce_sc3_magic,gene_filter = FALSE,
              ks=3:8,pct_dropout_min = 10, pct_dropout_max = 95,k_estimator=T,rand_seed=66666,biology=T)

pdf("pheatmap_magic_sc3.pdf",width = 10,height = 20)
sc3_plot_markers(sce_sc3_magic,k=8,show_pdata=c("bias_cDC1","bias_cDC2","bias_pDC"),fontsize=5,p.val=0.05,auroc=0.5)
dev.off()



#######

data = read.csv("/Users/tian.l/Dropbox/research/sis_seq/immgen/V1_ImmGenn-Official-Oct-2012.csv", skip=2, check.names=FALSE, as.is=TRUE)

data = data[data$`Gene Symbol` %in% rownames(sce_sc3_magic),]
data = data[!duplicated(data$`Gene Symbol`),]
rownames(data) = data$`Gene Symbol`
data = log2(data[,c("DC.8+.Sp.ST#6","DC.4+.Sp.ST#5","DC.pDC.8-.Sp#2","SC.CDP.BM#1")])
colnames(data) = c("cDC1_Immgen","cDC2_Immgen","pDC_Immgen","CDP_Immgen")

data = as.data.frame(t(scale(t(data))))

ann_colors = list(
  bias_cDC1 = c("#FFFFFF", "#FF7272"),
  bias_cDC2 = c("#FFFFFF", "#74FF72"),
  bias_pDC = c("#FFFFFF", "#7276FF"),
  cDC1_Immgen=c("blue", "red"),
  cDC2_Immgen=c("blue", "red"),
  pDC_Immgen=c("blue", "red"),
  CDP_Immgen=c("blue","red")
)

new_dist = dist(data)+0.5*dist(logcounts(sce_sc3_magic)[rownames(data),])
comb_hc_row = hclust(new_dist)


new_dist = dist(as.data.frame(colData(sce_sc3_magic)[,c("bias_cDC1", "bias_cDC2", "bias_pDC")]))+
  1.2*dist(t(logcounts(sce_sc3_magic)[rownames(data),]))
comb_hc_col = hclust(new_dist)

pdf("heatmap_magic_biggg.pdf",width = 10,height = 50)
pheatmap(logcounts(sce_sc3_magic)[rownames(data),], 
         annotation_col=as.data.frame(colData(sce_sc3_magic)[,c("bias_cDC1", "bias_cDC2", "bias_pDC")]),
         annotation_row=data,
         cluster_cols = comb_hc_col,
         #cluster_rows = comb_hc_row,
         #annotation_row=as.data.frame(rowData(sce_sc3_magic)[,c("cDC1_Immgen","cDC2_Immgen","pDC_Immgen","CDP_Immgen")]),
         annotation_colors=ann_colors, show_colnames = FALSE,
         border_color=NA,fontsize=5)
dev.off()

pdf("heatmap_magic_sisseq.pdf",width = 5,height = 4)
pheatmap(logcounts(sce_sc3_magic)[rownames(data),], 
         annotation_col=as.data.frame(colData(sce_sc3_magic)[,c("bias_cDC1", "bias_cDC2", "bias_pDC")]),
         #annotation_row=data,
         #cluster_cols = comb_hc_col,
         #cluster_rows = comb_hc_row,
         show_rownames = F,
         #annotation_row=as.data.frame(rowData(sce_sc3_magic)[,c("cDC1_Immgen","cDC2_Immgen","pDC_Immgen","CDP_Immgen")]),
         annotation_colors=ann_colors, show_colnames = FALSE,
         border_color=NA,fontsize=5)
dev.off()


new_dist = dist(data)+0.5*dist(log2(normcounts(sce_sc3_magic)[rownames(data),]+1))
comb_hc_row = hclust(new_dist)


new_dist = dist(as.data.frame(colData(sce_sc3_magic)[,c("bias_cDC1", "bias_cDC2", "bias_pDC")]))+
  1.5*dist(t(log2(normcounts(sce_sc3_magic)[rownames(data),]+1)))
comb_hc_col = hclust(new_dist)

pdf("heatmap_no_magic_biggg.pdf",width = 10,height = 50)
pheatmap(log2(normcounts(sce_sc3_magic)[rownames(data),]+1), 
         annotation_col=as.data.frame(colData(sce_sc3_magic)[,c("bias_cDC1", "bias_cDC2", "bias_pDC")]),
         annotation_row=data,
         cluster_cols = comb_hc_col,
         #cluster_rows = comb_hc_row,
         #annotation_row=as.data.frame(rowData(sce_sc3_magic)[,c("cDC1_Immgen","cDC2_Immgen","pDC_Immgen","CDP_Immgen")]),
         annotation_colors=ann_colors, show_colnames = FALSE,
         border_color=NA,fontsize=5)
dev.off()


