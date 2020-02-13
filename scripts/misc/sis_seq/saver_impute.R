library(SAVER)
library(pheatmap)
library(Seurat)

all_mk_genes = c("Siglech","Bcor","Wdr5","Rag2","Ccr9","Csf1","Klf4","Batf3","Irf4","Mpo","Tcf4","Zeb2","Gfi1","Ltbr","Runx2","Zbtb46","Runx1","Irf8","Id2","Mycl","Cbfa2t3")

setwd("~/Dropbox/research/sis_seq")
saver_output = saver(allcounts1$counts,ncores=1,
                     pred.genes=which(rownames(allcounts1$counts) %in% c(pDC_genes, cDC1_genes, cDC2_genes)),estimates.only=TRUE)
saveRDS(saver_output,file="saver_output.Rds")

saver_sel = log2(saver_output[rownames(saver_output) %in% c(pDC_genes, cDC1_genes, cDC2_genes),]+1)
saver_sel = t(scale(t(saver_sel)))
saver_sel[saver_sel<(-2.5)] = -2.5
saver_sel[saver_sel>2.5] = 2.5
ann_colors = list(
  bias_cDC1 = c("#FFFFFF", "#FF7272"),
  bias_cDC2 = c("#FFFFFF", "#74FF72"),
  bias_pDC = c("#FFFFFF", "#7276FF"),
  cDC1_Immgen=c("blue", "red"),
  cDC2_Immgen=c("blue", "red"),
  pDC_Immgen=c("blue", "red"),
  CDP_Immgen=c("blue","red")
)
mk_rowname = rownames(saver_sel)
mk_rowname[!(mk_rowname %in% all_mk_genes)] = ""
tp_fate = as.data.frame(fate_anno_pDC[,c("bias_cDC1", "bias_cDC2", "bias_pDC")])
rownames(tp_fate) = fate_anno_pDC$cell_name
pheatmap(saver_sel, 
         annotation_col=tp_fate,
         #annotation_row=data,
         #cluster_cols = comb_hc_col,
         #cluster_rows = comb_hc_row,
         #scale="column",
         labels_row=mk_rowname,
         show_rownames = T,
         color=PurpleAndYellow(),
         #annotation_row=as.data.frame(rowData(sce_sc3_magic)[,c("cDC1_Immgen","cDC2_Immgen","pDC_Immgen","CDP_Immgen")]),
         annotation_colors=ann_colors, show_colnames = FALSE,
         border_color=NA,fontsize=8,
         filename="SAVER_all_fate_gene_heatmap.pdf",
         width = 8,height = 7)
