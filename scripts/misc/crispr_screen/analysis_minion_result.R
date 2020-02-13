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


bc_list <-lapply(str_pad(1:12, 2, pad = "0"), function(x){read.csv(paste0("~/Dropbox/research/sis_seq/validation/CRISPR_screen/long_read/cnt/bc",x,".csv"), row.names=1)})

gRNA_lib <- read.csv("~/Dropbox/research/sis_seq/selected_gRNA_broadgpp.brie.library_457X4genes.csv", stringsAsFactors=FALSE)
selected_genes_3fate <- read.csv("~/Dropbox/research/sis_seq/selected_genes_3fate.csv", stringsAsFactors=FALSE)
selected_genes_3fate = selected_genes_3fate[!duplicated(selected_genes_3fate$gene_names),]
rownames(selected_genes_3fate) = selected_genes_3fate$gene_names

combined_df = Reduce(cbind,lapply(bc_list,function(x){x[gRNA_lib$sgRNA.Target.Sequence,]}))
combined_df[is.na(combined_df)] = 0

pdf("aaa.pdf",width = 15,height = 15)
pairs(log2(combined_df+1))
dev.off()

colnames(combined_df) = c("WT.cDC1.A","WT.cDC1.B","Cas9.cDC1.A","Cas9.cDC1.B",
                          "WT.cDC2.A","WT.cDC2.B","Cas9.cDC2.A","Cas9.cDC2.B",
                          "WT.pDC.A","WT.pDC.B","Cas9.pDC.A","Cas9.pDC.B")

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

pdf("bbb.pdf",width = 15,height = 15)
pairs(log2(cnt_df+1))
dev.off()



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
dt1 = decideTests(lrt1,lfc = 1,p.value = 0.1,adjust.method="fdr")
summary(dt1)

lrt2 = glmLRT(fit,contrast=contr.matrix[,2])
dt2 = decideTests(lrt2,lfc = 1,p.value = 0.1,adjust.method="fdr")
tp2 = topTags(lrt2,n=Inf,adjust.method="fdr",sort.by="none")
summary(dt2)

lrt3 = glmLRT(fit,contrast=contr.matrix[,3])
tp3 = topTags(lrt3,n=Inf,adjust.method="fdr",sort.by="none")
dt3 = decideTests(lrt3,lfc = 1,p.value = 0.1,adjust.method="fdr")
summary(dt3)


library(Glimma)

glMDPlot(lrt1, status=dt1, main=colnames(tfit)[1],anno=dge$genes, display.columns=c("gene_names","overlap","cDC1_FDR","cDC2_FDR","pDC_FDR"),counts=dge$counts, groups=K,html="cDC1 WT vs Cas9")
glMDPlot(lrt2, status=dt2, main=colnames(tfit)[2],anno=dge$genes, display.columns=c("gene_names","overlap","cDC1_FDR","cDC2_FDR","pDC_FDR"),counts=dge$counts , groups=K,html="cDC2 WT vs Cas9")
glMDPlot(lrt3, status=dt3, main=colnames(tfit)[3],anno=dge$genes, display.columns=c("gene_names","overlap","cDC1_FDR","cDC2_FDR","pDC_FDR"),counts=dge$counts, groups=K,html="pDC WT vs Cas9")


DE_FC = data.frame(gene_names=tp1@.Data[[1]]$gene_names,
                 cDC1_logFC=tp1@.Data[[1]]$logFC,
                 cDC1_FDR=tp1@.Data[[1]]$FDR,
                 cDC2_FDR=tp2@.Data[[1]]$FDR,
                 pDC_FDR=tp3@.Data[[1]]$FDR,
                 cDC2_logFC=tp2@.Data[[1]]$logFC,
                 pDC_logFC=tp3@.Data[[1]]$logFC)
DE_FC$min_FDR = as.numeric(apply(DE_FC,1,function(x){min(x[3:5])}))

library(data.table)
DE_FC = setDT(DE_FC)[ , .SD[which.min(min_FDR)], by = gene_names]
DE_FC = as.data.frame(DE_FC)
rownames(DE_FC) = DE_FC$gene_names

DE_FC = DE_FC[DE_FC$min_FDR<0.1,]

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

htmlwidgets::saveWidget(widget=aa,"sig_genes_logFC_ternary.html", selfcontained = TRUE)
#pty_search(aa,defaultSearchVal="Irf8", SepHTMLname="CRISPR_logFC.html")
library(tidyr)
library(ggpubr)

plot_gene_box = function(gene_name){
  g1 = all_df[all_df$Target.Gene.Symbol == gene_name,colnames(combined_df)]
  g1 <- gather(g1, group, read_number, 1:ncol(g1), factor_key=TRUE)
  g1$condition = "WT"
  g1$condition[grepl("Cas9",g1$group)] = "Cas9"
  g1$cell_type = "cDC1"
  g1$cell_type[grepl("cDC2",g1$group)] = "cDC2"
  g1$cell_type[grepl("pDC",g1$group)] = "pDC"
  
  
  p1 = ggplot(data=g1,aes(x=cell_type,y=read_number,fill=condition))+
    geom_boxplot()+
    stat_compare_means(aes(group = condition), label = "p.signif")+
    ggtitle(gene_name)+
    theme_bw()+
    theme(text = element_text(size=15))
  p1
}

pdf("CRISPR_MinION_all_boxplot.pdf",width = 6,height = 4)
a = c("Batf3","Jade1","Bcor","Id2","Tcf4","Spib","Klf4","Siglech","Ccl2","Rag2","Ccr9","Csf1","Clec11a","Mpo")
for (a_g in unique(c(a,rownames(DE_FC)))) {
  if (a_g %in% all_df$Target.Gene.Symbol){
    print(plot_gene_box(a_g))
  }
}
dev.off()


################


dd = read.csv(text="Lcn2
Mlf2
Tfec
Neb
Bcor
Mycl
Mycn
Ifi204
Nsg1
Jade1
Anxa1
Anxa2
Fn1
Zfp652
Pilra
Txk
Zbtb14
Pqlc3
Runx1
Runx2
Sh2d2a
Plac8
Cux2",header=F,stringsAsFactors=F)

def = dd$V1
def


dd = read.csv(text="Lyn
Pou2af1
Zfp263
Tyrobp
Zfhx2
Prrg2
Pacsin2
Pygl
Ctnnbip1
Havcr2",header=F,stringsAsFactors=F)

may = dd$V1
may



broadgpp.brie.library.contents <- read.delim("~/Dropbox/research/sis_seq/broadgpp-brie-library-contents.txt", stringsAsFactors=FALSE)

def[!(def %in% broadgpp.brie.library.contents$Target.Gene.Symbol)]

def_gRNA = broadgpp.brie.library.contents[broadgpp.brie.library.contents$Target.Gene.Symbol %in% def,]


res_seq = c()
for (i in def)
{
  if (any(dt1[dge$genes$gene_names == i] !=0)){
    res_seq = c(res_seq, all_df$Target.Context.Sequence[(dt1!=0) & (dge$genes$gene_names == i)])
  }
  else if (any(dt2[dge$genes$gene_names == i] !=0)){
    res_seq = c(res_seq, all_df$Target.Context.Sequence[(dt2!=0) & (dge$genes$gene_names == i)])
  }
  else if (any(dt3[dge$genes$gene_names == i] !=0)){
    res_seq = c(res_seq, all_df$Target.Context.Sequence[(dt3!=0) & (dge$genes$gene_names == i)])
  }
  else{
    print(i)
    res_seq = c(res_seq, all_df$Target.Context.Sequence[(dge$genes$gene_names == i)])
  }
}

def_gRNA_sel = def_gRNA[def_gRNA$Target.Context.Sequence %in% res_seq,]
def_list = list()
for (i in def){
  if (dim(def_gRNA_sel[def_gRNA_sel$Target.Gene.Symbol == i,])[1]>=2)
  {
    def_list[[i]] = (def_gRNA_sel[def_gRNA_sel$Target.Gene.Symbol == i,][1:2,])
  }
  else
  {
    tp = def_gRNA_sel[def_gRNA_sel$Target.Gene.Symbol == i,]
    tp1 = def_gRNA[def_gRNA$Target.Gene.Symbol==tp$Target.Gene.Symbol & def_gRNA$sgRNA.Target.Sequence != tp$sgRNA.Target.Sequence,][1,]
    def_list[[i]] = (rbind(tp,tp1))
  }
}

def_new = Reduce(rbind,def_list)

########

may[!(may %in% broadgpp.brie.library.contents$Target.Gene.Symbol)]

may_gRNA = broadgpp.brie.library.contents[broadgpp.brie.library.contents$Target.Gene.Symbol %in% c(may,"Batf3","Irf8","Tcf4"),]





res_seq = c()
for (i in c(may,"Batf3","Irf8","Tcf4"))
{
  if (any(dt1[dge$genes$gene_names == i] !=0)){
    res_seq = c(res_seq, all_df$Target.Context.Sequence[(dt1!=0) & (dge$genes$gene_names == i)])
  }
  else if (any(dt2[dge$genes$gene_names == i] !=0)){
    res_seq = c(res_seq, all_df$Target.Context.Sequence[(dt2!=0) & (dge$genes$gene_names == i)])
  }
  else if (any(dt3[dge$genes$gene_names == i] !=0)){
    res_seq = c(res_seq, all_df$Target.Context.Sequence[(dt3!=0) & (dge$genes$gene_names == i)])
  }
  else{
    print(i)
  }
}

may_gRNA = may_gRNA[may_gRNA$Target.Context.Sequence %in% res_seq,]
may_gRNA = may_gRNA[!duplicated(may_gRNA$Target.Gene.Symbol),]
may_gRNA

final_df = rbind(def_new, may_gRNA)

write.csv(final_df,file="selected_gRNA_SIS-seq_screen_v1.csv")


