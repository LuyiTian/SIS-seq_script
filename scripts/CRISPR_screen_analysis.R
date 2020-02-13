library(stringr)
pty_search = function (plotlyOut, defaultSearchVal="Bcor", SepHTMLname="plotly_with_search") 
{
  require(plotly)
  require(htmlwidgets)
  require(htmltools)
  
  p <- htmlwidgets::appendContent(plotlyOut, htmltools::tags$input(id='inputText', value=defaultSearchVal, ''), htmltools::tags$button(id='buttonSearch', 'Search'))
  p <- htmlwidgets::appendContent(p, htmltools::tags$script(HTML(
    'document.getElementById("buttonSearch").addEventListener("click", function()
  {        
  var i = 0;
  var j = 0;
  var found = [];
  var myDiv = document.getElementsByClassName("js-plotly-plot")[0]
  var data = JSON.parse(document.querySelectorAll("script[type=\'application/json\']")[0].innerHTML);
  for (i = 0 ;i < data.x.data.length; i += 1) {
  for (j = 0; j < data.x.data[i].text.length; j += 1) {
  if (data.x.data[i].text[j].indexOf(document.getElementById("inputText").value) !== -1) {
  found.push({curveNumber: i, pointNumber: j});
  }
  }
  }
  Plotly.Fx.hover(myDiv, found);
  }  
);'))) 
  htmlwidgets::saveWidget(p, SepHTMLname)
  p
}  

str_pad(1:12, 2, pad = "0")


bc_list <-lapply(str_pad(1:12, 2, pad = "0"), function(x){read.csv(paste0("~/Dropbox/research/sis_seq/validation/CRISPR_screen/long_read/cnt_p/bc",x,".csv"), row.names=1)})

gRNA_lib <- read.csv("~/Dropbox/research/sis_seq/selected_gRNA_broadgpp.brie.library_457X4genes.csv", stringsAsFactors=FALSE)
selected_genes_3fate <- read.csv("~/Dropbox/research/sis_seq/selected_genes_3fate.csv", stringsAsFactors=FALSE)
selected_genes_3fate = selected_genes_3fate[!duplicated(selected_genes_3fate$gene_names),]
rownames(selected_genes_3fate) = selected_genes_3fate$gene_names

combined_df = Reduce(cbind,lapply(bc_list,function(x){x[gRNA_lib$sgRNA.Target.Sequence,]}))
combined_df[is.na(combined_df)] = 0



colnames(combined_df) = c("WT.cDC1.A","WT.cDC1.B","Cas9.cDC1.A","Cas9.cDC1.B",
                          "WT.cDC2.A","WT.cDC2.B","Cas9.cDC2.A","Cas9.cDC2.B",
                          "WT.pDC.A","WT.pDC.B","Cas9.pDC.A","Cas9.pDC.B")

plot(colSums(combined_df))

all_df = cbind(gRNA_lib, combined_df)
all_df = all_df[order(all_df$Target.Gene.Symbol),]

selected_genes = selected_genes_3fate[all_df$Target.Gene.Symbol,]
selected_genes$gene_names = all_df$Target.Gene.Symbol
selected_genes = selected_genes[order(selected_genes$gene_names),]
cnt_df = all_df[,colnames(combined_df)]

all_df[all_df$Target.Gene.Symbol=="Pilra",] 

tmp_df = merge(selected_genes_3fate,all_df,by.x = "gene_names", by.y="Target.Gene.Symbol",all.x = FALSE,all.y=TRUE)
tmp_df = tmp_df[order(tmp_df$gene_names),]
#cnt_gene_df = aggregate(. ~ Target.Gene.Symbol, data = cnt_df, sum)

#rownames(cnt_df) = all_df$Target.Gene.Symbol
#cnt_gene_df = cnt_gene_df[,-1]





library(edgeR)
library(limma)

dge <- DGEList(counts=cnt_df,genes=selected_genes)
dge <- calcNormFactors(dge)

conditions = rep(c(rep("WT",2),rep("cas9",2)),3)
DC_type = rep(c("cDC1","cDC2","pDC"),each=4)

K = paste(conditions,DC_type,sep="_")

design_mat = model.matrix(~ 0+K)


contr.matrix <- makeContrasts(
  cDC1 = Kcas9_cDC1-KWT_cDC1,
  cDC2 = Kcas9_cDC2 - KWT_cDC2, 
  pDC = Kcas9_pDC - KWT_pDC,
  levels = colnames(design_mat))
contr.matrix


# v <- voom(dge, design_mat, plot=TRUE)
# 
# 
# vfit <- lmFit(v, design_mat)
# vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
# efit <- eBayes(vfit)
# plotSA(efit)
# 
# tfit <- treat(vfit, lfc=1)
# dt <- decideTests(tfit,lfc = 1,p.value = 1)
# summary(dt)

all_cnt = estimateDisp(dge, design_mat, robust=TRUE)

fit = glmFit(all_cnt, design_mat)
# Run likelihood ratio tests
lrt1 = glmLRT(fit,contrast=contr.matrix[,1])
tp1 = topTags(lrt1,n=Inf,adjust.method="fdr",sort.by="none")
dt1 = decideTests(lrt1,lfc = 0.5,p.value = 0.01,adjust.method="fdr")
summary(dt1)

lrt2 = glmLRT(fit,contrast=contr.matrix[,2])
dt2 = decideTests(lrt2,lfc = 0.5,p.value = 0.01,adjust.method="fdr")
tp2 = topTags(lrt2,n=Inf,adjust.method="fdr",sort.by="none")
summary(dt2)

lrt3 = glmLRT(fit,contrast=contr.matrix[,3])
tp3 = topTags(lrt3,n=Inf,adjust.method="fdr",sort.by="none")
dt3 = decideTests(lrt3,lfc = 0.5,p.value = 0.01,adjust.method="fdr")
summary(dt3)

#### roast gene set test
idx_list = list()
for (i in 1:length(unique(dge$genes$gene_names)))
{
  idx_list[[unique(dge$genes$gene_names)[i]]] = (4*i-3):(4*i)
}

r_res1 = mroast(all_cnt, idx_list, design_mat, contrast=contr.matrix[,1])
r_res2 = mroast(all_cnt, idx_list, design_mat, contrast=contr.matrix[,2])
r_res3 = mroast(all_cnt, idx_list, design_mat, contrast=contr.matrix[,3])


r_res1_down = r_res1[r_res1$PropDown>=0.5 & r_res1$PropUp==0,]
r_res2_down = r_res2[r_res2$PropDown>=0.5 & r_res2$PropUp==0,]
r_res3_down = r_res3[r_res3$PropDown>=0.5 & r_res3$PropUp==0,]

res1_d = setdiff(rownames(r_res1_down), union(rownames(r_res2_down),rownames(r_res3_down)))
r_res1_down_only = r_res1_down[res1_d,]
r_res1_down_only = r_res1_down_only[rownames(r_res1_down_only) %in% dge$genes$gene_names[dt1==-1],]
r_res1_down_only

write.csv(r_res1_down_only,file="only_cDC1_down_roast.csv")

res2_d = setdiff(rownames(r_res2_down), union(rownames(r_res1_down),rownames(r_res3_down)))
r_res2_down_only = r_res2_down[res2_d,]
r_res2_down_only = r_res2_down_only[rownames(r_res2_down_only) %in% dge$genes$gene_names[dt2==-1],]
r_res2_down_only

write.csv(r_res2_down_only,file="only_cDC2_down_roast.csv")

res3_d = setdiff(rownames(r_res3_down), union(rownames(r_res1_down),rownames(r_res2_down)))
r_res3_down_only = r_res3_down[res3_d,]
r_res3_down_only = r_res3_down_only[rownames(r_res3_down_only) %in% dge$genes$gene_names[dt3==-1],]
r_res3_down_only

write.csv(r_res3_down_only,file="only_pDC_down_roast.csv")
####


DE_FC = data.frame(gene_names=tp1@.Data[[1]]$gene_names,
                 cDC1_logFC=tp1@.Data[[1]]$logFC,
                 cDC1_FDR=tp1@.Data[[1]]$FDR,
                 cDC2_FDR=tp2@.Data[[1]]$FDR,
                 pDC_FDR=tp3@.Data[[1]]$FDR,
                 cDC2_logFC=tp2@.Data[[1]]$logFC,
                 pDC_logFC=tp3@.Data[[1]]$logFC)
DE_FC$min_FDR = as.numeric(apply(DE_FC,1,function(x){min(x[3:5])}))

DE_FC_sig = DE_FC[DE_FC$min_FDR<0.01,]
row_anno = data.frame(down = ifelse(DE_FC_sig$gene_names %in% c(rownames(r_res1_down_only),
                                                        rownames(r_res2_down_only),
                                                        rownames(r_res3_down_only)),"down","others"))
rownames(row_anno) = rownames(DE_FC_sig)
#row_anno$down = as.factor(row_anno$down)
pdf("heatmap_sig_gRNA_LogFoldChange.pdf",height = 40,width = 10)
pheatmap::pheatmap(DE_FC_sig[,c("cDC1_logFC","cDC2_logFC","pDC_logFC")],
                   annotation_row=row_anno,
                   annotation_colors=list(down = c(down = "firebrick", others = "grey")),
                   scale = "column",
                   fontsize = 4,
                   labels_row=DE_FC_sig$gene_names)
dev.off()

all_down = (DE_FC$cDC1_logFC<0 & DE_FC$cDC2_logFC<0 & DE_FC$pDC_logFC <0) & (DE_FC$cDC1_FDR< 0.01 & DE_FC$cDC2_FDR< 0.01 & DE_FC$pDC_FDR < 0.01)
all_up = (DE_FC$cDC1_logFC>0 & DE_FC$cDC2_logFC>0 & DE_FC$pDC_logFC >0) & (DE_FC$cDC1_FDR< 0.01 & DE_FC$cDC2_FDR< 0.01 & DE_FC$pDC_FDR < 0.01)

cDC1_only = (DE_FC$cDC1_FDR< 0.01 & DE_FC$cDC2_FDR> 0.05 & DE_FC$pDC_FDR > 0.05)
cDC2_only = (DE_FC$cDC2_FDR< 0.01 & DE_FC$cDC1_FDR> 0.05 & DE_FC$pDC_FDR > 0.05)
pDC_only = (DE_FC$pDC_FDR< 0.01 & DE_FC$cDC2_FDR> 0.05 & DE_FC$cDC1_FDR > 0.05)



# ####
# library(RISmed)
# search_str1 = "((_GE_[Title/Abstract]) AND dendritic cells[Title/Abstract]) AND (development[Title/Abstract] OR differentiation[Title/Abstract])"
# search_str2 = "_GE_[Title/Abstract]"
# ff1 =function(x){QueryCount(EUtilsSummary(gsub("_GE_",x,search_str1)))}
# ff2 =function(x){QueryCount(EUtilsSummary(gsub("_GE_",x,search_str2)))}
# 
# ref_cnt1 = c()
# ref_cnt2 = c()
# for (g in unique(dge$genes$gene_names)){
#   print(g)
#   ref_cnt1 = c(ref_cnt1, ff1(g))
#   ref_cnt2 = c(ref_cnt2, ff2(g))
#   Sys.sleep(5)
# }
# names(ref_cnt1) = unique(dge$genes$gene_names)
# names(ref_cnt2) = unique(dge$genes$gene_names)
# 
# ####

###
library(Glimma)
K_l = factor(c(K,"NA"), levels = c("WT_cDC1", "cas9_cDC1", "WT_cDC2","cas9_cDC2","WT_pDC","cas9_pDC","NA"))
mtx = edgeR::cpm(dge$counts)
mtx = cbind(mtx,rep(0,nrow(mtx)))

gene_anno = dge$genes
gene_anno$specificity = NA

gene_anno$specificity[all_down] = "all_down"
gene_anno$specificity[all_up] = "all_up"
gene_anno$specificity[cDC1_only] = "cDC1_only"
gene_anno$specificity[cDC2_only] = "cDC2_only"
gene_anno$specificity[pDC_only] = "pDC_only"

glMDPlot(lrt1, status=dt1, main="cDC1",anno=gene_anno, display.columns=c("gene_names","overlap","regulatory_genes","specificity"),counts=mtx, transform=F,groups=K_l,html="p cDC1 WT vs Cas9")
glMDPlot(lrt2, status=dt2, main="cDC2",anno=gene_anno, display.columns=c("gene_names","overlap","regulatory_genes","specificity"),counts=mtx, transform=F,groups=K_l,html="p cDC2 WT vs Cas9")
glMDPlot(lrt3, status=dt3, main="pDC",anno=gene_anno, display.columns=c("gene_names","overlap","regulatory_genes","specificity"),counts=mtx, transform=F,groups=K_l,html="p pDC WT vs Cas9")

###

library(data.table)
DE_FC = setDT(DE_FC)[ , .SD[which.min(min_FDR)], by = gene_names]
DE_FC = as.data.frame(DE_FC)
rownames(DE_FC) = DE_FC$gene_names

DE_FC = DE_FC[DE_FC$min_FDR<0.001,]

pdf("heatmap_gene_LogFoldChange.pdf",height = 23,width = 5)
pheatmap::pheatmap(DE_FC[,c("cDC1_logFC","cDC2_logFC","pDC_logFC")],scale = "column",fontsize = 5)
dev.off()

tmp = DE_FC[,c("cDC1_logFC","cDC2_logFC","pDC_logFC")]
pairs(tmp)


range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#tmp = range01(sign(scale(tmp))*abs(scale(tmp))^(1/4))
tmp = as.data.frame(tmp/rowSums(tmp))

col <- grDevices::rgb(tmp$cDC1_logFC, tmp$cDC2_logFC, tmp$pDC_logFC,alpha=0.6)

tmp_sel = tmp[rownames(tmp) %in% c("Bcor","Zbtb46","Naaa","Tcf4","Irf8","Ccr9"),]

library(Ternary)
pdf("CRISPR_genes_logFC_ternary.pdf")
TernaryPlot('pDC_logFC', 'cDC2_logFC', 'cDC1_logFC', 
            grid.lines=5, grid.lty='dotted',
            grid.minor.lines = 1, grid.minor.lty='dotted')
AddToTernary(points,tmp[,c("pDC_logFC","cDC2_logFC","cDC1_logFC")],pch=20, cex=1,col=col)
AddToTernary(text, tmp_sel[,c("pDC_logFC","cDC2_logFC","cDC1_logFC")], rownames(tmp_sel), cex=0.8, font=2)
dev.off()


library(plotly)
tmp$gene_name = rownames(tmp)
aa = plot_ly(a=~pDC_logFC,
             b=~cDC1_logFC,
             c=~cDC2_logFC,
             text=~gene_name,
             data=tmp,
             type = 'scatterternary',
             mode = 'markers',
             marker = list(size = 5))%>%
  layout(
    title = "significant genes logFC ternary plot",
    ternary = list(
      aaxis = list(title = "pDC"),
      baxis = list(title = "cDC1"),
      caxis = list(title = "cDC2")
    ))
aa

htmlwidgets::saveWidget(widget=aa,"p_sig_genes_logFC_ternary.html", selfcontained = TRUE)
#pty_search(aa,defaultSearchVal="Irf8", SepHTMLname="CRISPR_logFC.html")
library(tidyr)
library(ggpubr)

plot_gene_box = function(gene_name){
  g1 = as.data.frame(edgeR::cpm(dge$counts)[dge$genes$gene_names == gene_name,])
  g1 <- gather(g1, group, read_number, 1:ncol(g1), factor_key=TRUE)
  g1$condition = "WT"
  g1$condition[grepl("Cas9",g1$group)] = "Cas9"
  g1$cell_type = "cDC1"
  g1$cell_type[grepl("cDC2",g1$group)] = "cDC2"
  g1$cell_type[grepl("pDC",g1$group)] = "pDC"
  
  
  p1 = ggplot(data=g1,aes(x=cell_type,y=read_number,fill=condition))+
    geom_boxplot(coef = 500)+
    geom_point(size=1, position = position_jitterdodge())+
    stat_compare_means(aes(group = condition), label = "p.signif")+
    ggtitle(gene_name)+
    theme_bw()+
    theme(text = element_text(size=15))
  p1
}

pdf("CRISPR_Promethion_all_boxplot.pdf",width = 6,height = 4)
a = c("Batf3","Jade1","Bcor","Id2","Tcf4","Spib","Klf4","Siglech","Ccl2","Rag2","Ccr9","Csf1","Clec11a","Mpo")
for (a_g in unique(c(a,rownames(r_res1_down_only),rownames(r_res2_down_only),rownames(r_res3_down_only)))) {
  if (a_g %in% all_df$Target.Gene.Symbol){
    print(plot_gene_box(a_g))
  }
}
dev.off()


