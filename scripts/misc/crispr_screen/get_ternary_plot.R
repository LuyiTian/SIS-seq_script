

prop_WT = as.data.frame(mtx[,grepl("WT",colnames(mtx))])

prop_WT$cDC1_WT = prop_WT$WT.cDC1.A+prop_WT$WT.cDC1.B
prop_WT$cDC2_WT = prop_WT$WT.cDC2.A+prop_WT$WT.cDC2.B
prop_WT$pDC_WT = prop_WT$WT.pDC.A+prop_WT$WT.pDC.B

prop_WT$sum = rowSums(prop_WT[,c("cDC1_WT","cDC2_WT","pDC_WT")])

prop_WT$cDC1_WT = prop_WT$cDC1_WT/prop_WT$sum
prop_WT$cDC2_WT = prop_WT$cDC2_WT/prop_WT$sum
prop_WT$pDC_WT = prop_WT$pDC_WT/prop_WT$sum

library(ggplot2)
library(ggtern)

pdf("WT_VS_Cas9.pdf")
ggtern(data=prop_WT,aes(cDC1_WT,cDC2_WT,pDC_WT)) + 
  geom_point(alpha=0.2,size=0.5)+
  stat_density_tern()+
  theme_bw()

WT_sig_val = list()
WT_sig_val$cDC1_lo = unname(quantile(prop_WT$cDC1_WT,0.05))
WT_sig_val$cDC1_hi = unname(quantile(prop_WT$cDC1_WT,0.95))
WT_sig_val$cDC2_lo = unname(quantile(prop_WT$cDC2_WT,0.05))
WT_sig_val$cDC2_hi = unname(quantile(prop_WT$cDC2_WT,0.95))
WT_sig_val$pDC_lo =  unname(quantile(prop_WT$pDC_WT,0.05))
WT_sig_val$pDC_hi =  unname(quantile(prop_WT$pDC_WT,0.95))



prop_Cas9 = as.data.frame(mtx[,grepl("Cas9",colnames(mtx))])

prop_Cas9$cDC1_Cas9 = prop_Cas9$Cas9.cDC1.A+prop_Cas9$Cas9.cDC1.B
prop_Cas9$cDC2_Cas9 = prop_Cas9$Cas9.cDC2.A+prop_Cas9$Cas9.cDC2.B
prop_Cas9$pDC_Cas9 = prop_Cas9$Cas9.pDC.A+prop_Cas9$Cas9.pDC.B

prop_Cas9$sum = rowSums(prop_Cas9[,c("cDC1_Cas9","cDC2_Cas9","pDC_Cas9")])

prop_Cas9$cDC1_Cas9 = prop_Cas9$cDC1_Cas9/prop_Cas9$sum
prop_Cas9$cDC2_Cas9 = prop_Cas9$cDC2_Cas9/prop_Cas9$sum
prop_Cas9$pDC_Cas9 = prop_Cas9$pDC_Cas9/prop_Cas9$sum

keep_cpm = prop_Cas9$sum>36
prop_Cas9 = prop_Cas9[keep_cpm,] # avg cpm > 6

ggtern(data=prop_Cas9,aes(cDC1_Cas9,cDC2_Cas9,pDC_Cas9)) + 
  geom_point(alpha=0.2,size=0.5)+
  stat_density_tern()+
  theme_bw()
dev.off()


#i  = unique(dge$genes$gene_names)[1]
gRNA_fate_df = NULL
for (i in unique(dge$genes$gene_names)){
  tmp_v = rep(FALSE,6) 
  prop_mtx=prop_Cas9[dge$genes$gene_name == i,c("cDC1_Cas9","cDC2_Cas9","pDC_Cas9")]
  cDC1_lo = unname(table(prop_mtx[,1] < WT_sig_val$cDC1_lo)["TRUE"])
  cDC1_hi = unname(table(prop_mtx[,1] > WT_sig_val$cDC1_hi)["TRUE"])

  if (is.na(cDC1_hi) & !is.na(cDC1_lo)){
    if (cDC1_lo>1){tmp_v[1] = TRUE}
  }else if(!is.na(cDC1_hi) & is.na(cDC1_lo)){
    if (cDC1_hi>1){tmp_v[2] = TRUE}
  }
  
  cDC2_lo = unname(table(prop_mtx[,2] < WT_sig_val$cDC2_lo)["TRUE"])
  cDC2_hi = unname(table(prop_mtx[,2] > WT_sig_val$cDC2_hi)["TRUE"])
  if (is.na(cDC2_hi) & !is.na(cDC2_lo)){
    if (cDC2_lo>1){tmp_v[3] = TRUE}
  }else if(!is.na(cDC2_hi) & is.na(cDC2_lo)){
    if (cDC2_hi>1){tmp_v[4] = TRUE}
  }
  
  pDC_lo = unname(table(prop_mtx[,3] < WT_sig_val$pDC_lo)["TRUE"])
  pDC_hi = unname(table(prop_mtx[,3] > WT_sig_val$pDC_hi)["TRUE"])
  if (is.na(pDC_hi) & !is.na(pDC_lo)){
    if (pDC_lo>1){tmp_v[5] = TRUE}
  }else if(!is.na(pDC_hi) & is.na(pDC_lo)){
    if (pDC_hi>1){tmp_v[6] = TRUE}
  }
  
  if(is.null(gRNA_fate_df)){
    gRNA_fate_df = data.frame(gene_name=i,cDC1_lo=tmp_v[1],cDC1_hi=tmp_v[2],
                              cDC2_lo=tmp_v[3],cDC2_hi=tmp_v[4],
                              pDC_lo=tmp_v[5],pDC_hi=tmp_v[6])
  }else{
    gRNA_fate_df = rbind(gRNA_fate_df,data.frame(gene_name=i,cDC1_lo=tmp_v[1],cDC1_hi=tmp_v[2],
                                                 cDC2_lo=tmp_v[3],cDC2_hi=tmp_v[4],
                                                 pDC_lo=tmp_v[5],pDC_hi=tmp_v[6]))
  }
}

for(i in 2:7){gRNA_fate_df[,i] = as.logical(gRNA_fate_df[,i])}

gRNA_fate_df$keep = apply(gRNA_fate_df,1,function(x){table(x[2:7])["FALSE"]<6})

gRNA_fate_df = gRNA_fate_df[gRNA_fate_df$keep,]

prop_Cas9_sel = prop_Cas9
prop_Cas9_sel$gene_name = dge$genes$gene_name[keep_cpm]
write.csv(prop_Cas9_sel,file="gRNA_fate_full_table.csv")
prop_Cas9_sel = prop_Cas9_sel[prop_Cas9_sel$gene_name %in% gRNA_fate_df$gene_name,]

keep1 = prop_Cas9_sel$cDC1_Cas9 < WT_sig_val$cDC1_lo | prop_Cas9_sel$cDC1_Cas9 > WT_sig_val$cDC1_hi
keep2 = prop_Cas9_sel$cDC2_Cas9 < WT_sig_val$cDC2_lo | prop_Cas9_sel$cDC2_Cas9 > WT_sig_val$cDC2_hi
keep3 = prop_Cas9_sel$pDC_Cas9 < WT_sig_val$pDC_lo | prop_Cas9_sel$pDC_Cas9 > WT_sig_val$pDC_hi
prop_Cas9_sel = prop_Cas9_sel[(keep1|keep2)|keep3,]
library(plotly)
aa = plot_ly(a=~pDC_Cas9,
             b=~cDC1_Cas9,
             c=~cDC2_Cas9,
             text=~gene_name,
             data=prop_Cas9_sel,
             type = 'scatterternary',
             mode = 'markers',
             marker = list(size = 5))%>%
  layout(
    title = "significant gRNAs fate outcome",
    ternary = list(
      aaxis = list(title = "pDC"),
      baxis = list(title = "cDC1"),
      caxis = list(title = "cDC2")
    ))
aa
htmlwidgets::saveWidget(widget=aa,"sig_gRNA_fate_proportion.html", selfcontained = TRUE)

rownames(gRNA_fate_df) = gRNA_fate_df$gene_name
gRNA_fate_df = gRNA_fate_df[,-c(1,8)]
tmp = apply(gRNA_fate_df,2,as.numeric)
rownames(tmp) = rownames(gRNA_fate_df)
pdf("gene_fate_as_binary.pdf",width = 5,height = 10)
pheatmap::pheatmap(tmp,cluster_cols=FALSE)
dev.off()


dt1_new = dt1
dt1_new@.Data[,1] = 0
ix = rownames(prop_Cas9_sel)[prop_Cas9_sel$gene_name %in% rownames(gRNA_fate_df)[gRNA_fate_df$cDC1_lo]]
dt1_new[ix,1] = -1
ix = rownames(prop_Cas9_sel)[prop_Cas9_sel$gene_name %in% rownames(gRNA_fate_df)[gRNA_fate_df$cDC1_hi]]
dt1_new[ix,1] = 1

dt2_new = dt2
dt2_new@.Data[,1] = 0
ix = rownames(prop_Cas9_sel)[prop_Cas9_sel$gene_name %in% rownames(gRNA_fate_df)[gRNA_fate_df$cDC2_lo]]
dt2_new[ix,1] = -1
ix = rownames(prop_Cas9_sel)[prop_Cas9_sel$gene_name %in% rownames(gRNA_fate_df)[gRNA_fate_df$cDC2_hi]]
dt2_new[ix,1] = 1

dt3_new = dt3
dt3_new@.Data[,1] = 0
ix = rownames(prop_Cas9_sel)[prop_Cas9_sel$gene_name %in% rownames(gRNA_fate_df)[gRNA_fate_df$pDC_lo]]
dt3_new[ix,1] = -1
ix = rownames(prop_Cas9_sel)[prop_Cas9_sel$gene_name %in% rownames(gRNA_fate_df)[gRNA_fate_df$pDC_hi]]
dt3_new[ix,1] = 1

gene_anno$gRNA_id = rownames(gene_anno)

library(Glimma)
glMDPlot(lrt1, status=dt1_new, main="cDC1",anno=gene_anno, display.columns=c("gene_names","overlap","regulatory_genes","specificity","gRNA_id"),counts=mtx, transform=F,groups=K_l,html="prop cDC1 WT vs Cas9")
glMDPlot(lrt2, status=dt2_new, main="cDC2",anno=gene_anno, display.columns=c("gene_names","overlap","regulatory_genes","specificity","gRNA_id"),counts=mtx, transform=F,groups=K_l,html="prop cDC2 WT vs Cas9")
glMDPlot(lrt3, status=dt3_new, main="pDC",anno=gene_anno, display.columns=c("gene_names","overlap","regulatory_genes","specificity","gRNA_id"),counts=mtx, transform=F,groups=K_l,html="prop pDC WT vs Cas9")


library(ggpubr)
cond = c("Cas9.cDC1", "Cas9.cDC1", "Cas9.cDC2", "Cas9.cDC2", "Cas9.pDC", "Cas9.pDC",
         "WT.cDC1", "WT.cDC1", "WT.cDC2", "WT.cDC2", "WT.pDC", "WT.pDC")
pdf("sel_genes_fate_outcome.pdf",width = 10,height = 5)

for (i in unique(prop_Cas9_sel$gene_name)){
  ix = 1
  pp=list()
  tmp = prop_Cas9_sel[prop_Cas9_sel$gene_name==i,]
  ppp = ggtern(data=tmp,aes(pDC_Cas9,cDC1_Cas9,cDC2_Cas9)) + 
    geom_point(alpha=0.8,size=1)+
    stat_density_tern(data=prop_WT,aes(pDC_WT,cDC1_WT,cDC2_WT))+
    #Tlab("pDC") + Llab("cDC1") + Rlab("cDC2")+
    ggtitle(i)
  print(ppp)
  #pp[[ix]] = ggtern(data=tmp,aes(cDC1_Cas9,cDC2_Cas9,pDC_Cas9)) + 
  #  geom_point(alpha=0.8,size=1)+
  #  stat_density_tern(data=prop_WT,aes(cDC1_WT,cDC2_WT,pDC_WT))
  #ix = ix+1
  for (t in 1:nrow(tmp)){
    tmp_df = data.frame(condition=cond,norm_counts=unlist(c((tmp[t,1:6]),(prop_WT[rownames(tmp[t,]),1:6]))))
    tmp_df$condition = factor(tmp_df$condition, levels = c("WT.cDC1", "Cas9.cDC1", "WT.cDC2","Cas9.cDC2","WT.pDC","Cas9.pDC"))
    pp[[ix]] = ggplot(data=tmp_df,aes(x=condition,y=norm_counts))+geom_point()+theme_bw()+ggtitle(rownames(tmp[t,]))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ix = ix+1
  }
  print(ggarrange(plotlist=pp,nrow=1,ncol=length(pp)))
}
dev.off()



sel_genes = c("Zfp455","Ifi205","Zscan29","Tyrobp","Klf12","Ankrd24","Notch2","Nrxn1",
              "Zc3h12a","Adam19","Naaa","Trpa1","Celsr3","Mbd5","Ubr3","Drd1","Impdh2","Zfp652","Zfp831")

sel_snd = c("Havcr2","F2rl1","Il12rb1","Csf1r","Eif4ebp3","Tmem108","Wnt2b","Plxna2","Siglecg","Tespa1","Cfc1","Llgl1","Ccdc150")


library(ggpubr)
cond = c("Cas9.cDC1", "Cas9.cDC1", "Cas9.cDC2", "Cas9.cDC2", "Cas9.pDC", "Cas9.pDC",
         "WT.cDC1", "WT.cDC1", "WT.cDC2", "WT.cDC2", "WT.pDC", "WT.pDC")
pdf("sel_genes_top.pdf",width = 10,height = 5)

prop_Cas9_sel = prop_Cas9
prop_Cas9_sel$gene_name = dge$genes$gene_name[keep_cpm]

for (i in sel_genes){
  ix = 1
  pp=list()
  tmp = prop_Cas9_sel[prop_Cas9_sel$gene_name==i,]
  ppp = ggtern(data=tmp,aes(pDC_Cas9,cDC1_Cas9,cDC2_Cas9)) + 
    geom_point(alpha=0.8,size=1)+
    stat_density_tern(data=prop_WT,aes(pDC_WT,cDC1_WT,cDC2_WT))+
    #Tlab("pDC") + Llab("cDC1") + Rlab("cDC2")+
    ggtitle(i)
  print(ppp)
  #pp[[ix]] = ggtern(data=tmp,aes(cDC1_Cas9,cDC2_Cas9,pDC_Cas9)) + 
  #  geom_point(alpha=0.8,size=1)+
  #  stat_density_tern(data=prop_WT,aes(cDC1_WT,cDC2_WT,pDC_WT))
  #ix = ix+1
  for (t in 1:nrow(tmp)){
    tmp_df = data.frame(condition=cond,norm_counts=unlist(c((tmp[t,1:6]),(prop_WT[rownames(tmp[t,]),1:6]))))
    tmp_df$condition = factor(tmp_df$condition, levels = c("WT.cDC1", "Cas9.cDC1", "WT.cDC2","Cas9.cDC2","WT.pDC","Cas9.pDC"))
    pp[[ix]] = ggplot(data=tmp_df,aes(x=condition,y=norm_counts))+geom_point()+theme_bw()+ggtitle(rownames(tmp[t,]))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ix = ix+1
  }
  print(ggarrange(plotlist=pp,nrow=1,ncol=length(pp)))
}
dev.off()


pdf("sel_genes_second.pdf",width = 10,height = 5)

prop_Cas9_sel = prop_Cas9
prop_Cas9_sel$gene_name = dge$genes$gene_name[keep_cpm]

for (i in sel_snd){
  ix = 1
  pp=list()
  tmp = prop_Cas9_sel[prop_Cas9_sel$gene_name==i,]
  ppp = ggtern(data=tmp,aes(pDC_Cas9,cDC1_Cas9,cDC2_Cas9)) + 
    geom_point(alpha=0.8,size=1)+
    stat_density_tern(data=prop_WT,aes(pDC_WT,cDC1_WT,cDC2_WT))+
    #Tlab("pDC") + Llab("cDC1") + Rlab("cDC2")+
    ggtitle(i)
  print(ppp)
  #pp[[ix]] = ggtern(data=tmp,aes(cDC1_Cas9,cDC2_Cas9,pDC_Cas9)) + 
  #  geom_point(alpha=0.8,size=1)+
  #  stat_density_tern(data=prop_WT,aes(cDC1_WT,cDC2_WT,pDC_WT))
  #ix = ix+1
  for (t in 1:nrow(tmp)){
    tmp_df = data.frame(condition=cond,norm_counts=unlist(c((tmp[t,1:6]),(prop_WT[rownames(tmp[t,]),1:6]))))
    tmp_df$condition = factor(tmp_df$condition, levels = c("WT.cDC1", "Cas9.cDC1", "WT.cDC2","Cas9.cDC2","WT.pDC","Cas9.pDC"))
    pp[[ix]] = ggplot(data=tmp_df,aes(x=condition,y=norm_counts))+geom_point()+theme_bw()+ggtitle(rownames(tmp[t,]))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ix = ix+1
  }
  print(ggarrange(plotlist=pp,nrow=1,ncol=length(pp)))
}
dev.off()
