library(plotly)
library(Ternary)
pairs(all_list1$batch_DEs[,c("cDC1_logFC","cDC2_logFC","pDC_logFC")])


DE_FC = all_list1$batch_DEs

DE_FC = DE_FC[DE_FC$gene_names %in% c(cDC1_genes,cDC2_genes,pDC_genes),]
rownames(DE_FC) = DE_FC$gene_names

DE_FC = DE_FC[(DE_FC$cDC1_FDR<0.05 | DE_FC$cDC2_FDR<0.05) | DE_FC$pDC_FDR<0.05,]

plot_ly(x=DE_FC[,c("cDC1_logFC")],
        y=DE_FC[,c("cDC2_logFC")],
        z=DE_FC[,c("pDC_logFC")],
        sizes=1)


tmp = DE_FC[,c("cDC1_logFC","cDC2_logFC","pDC_logFC")]
colo = rep("others",length(rownames(tmp)))
colo[rownames(tmp)=="Bcor"] = "Bcor"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
tmp = range01(sign(scale(tmp))*abs(scale(tmp))^(1/3))
tmp = as.data.frame(tmp/rowSums(tmp))
aa = plot_ly(x=tmp[,c("cDC1_logFC")],
        y=tmp[,c("cDC2_logFC")],
        z=tmp[,c("pDC_logFC")],
        text=rownames(tmp),
        color=colo,
        colors=c("red","grey"),
        marker = list(size = 3))%>%
  layout(
    title = "significant genes logFC ternary plot",
    scene = list(
      xaxis = list(title = "cDC1"),
      yaxis = list(title = "cDC2"),
      zaxis = list(title = "pDC")
    ))
aa

htmlwidgets::saveWidget(widget=aa,"sig_genes_logFC_ternary.html", selfcontained = TRUE)

col <- grDevices::rgb(tmp$cDC1_logFC, tmp$cDC2_logFC, tmp$pDC_logFC,alpha=0.6)

tmp_sel = tmp[rownames(tmp) %in% c("Wdr5","Cbfa2t3","Bcor","Mycl","Batf3","Id2","Irf8","Runx1","Zbtb46","Tcf4","Spib","Runx2","Klf4","Ltbr","Siglech","Ccl2","Rag2","Gfi1","Zeb2","Ccr9","Tcf4","Csf1","Clec11a","Mpo"),]

pdf("sig_genes_logFC_ternary_big.pdf",width = 15,height = 15)
TernaryPlot('pDC_logFC', 'cDC2_logFC', 'cDC1_logFC', 
            grid.lines=5, grid.lty='dotted',
            grid.minor.lines = 1, grid.minor.lty='dotted')
AddToTernary(points,tmp[,c("pDC_logFC","cDC2_logFC","cDC1_logFC")],pch=20, cex=1,col=col)
AddToTernary(text, tmp_sel[,c("pDC_logFC","cDC2_logFC","cDC1_logFC")], rownames(tmp_sel), cex=0.8, font=2)
dev.off()
#library(ggtern)
#ggtern(data=tmp,aes(x=cDC1_logFC,y=cDC2_logFC,z=pDC_logFC))+geom_point()
