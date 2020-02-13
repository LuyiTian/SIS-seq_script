


tno_cnt = logcounts(sce)[,fate_anno_pDC$cell_name]

the_theme =   theme(text = element_text(size=15,face = "bold"),
                    plot.title = element_text(size=20,hjust = 0.5),
                    #axis.ticks.x=element_blank(),
                    #axis.text.x = element_blank(),
                    #axis.ticks.y=element_blank(),
                    #axis.text.y = element_blank(),
                    #panel.grid.major = element_blank(), 
                    #panel.grid.minor = element_blank(),
                    #panel.background = element_blank(),
                    axis.line = element_line(colour = "black"))

ggplot_cDC1_genes = names(glmout_cDC1)


p_list = list()
i=1

sel_mk_genes = c("Cbfa2t3","Mycl2","Id2","Batf3","Irf8","Runx1","Zbtb46","Runx2","Tcf4","Zeb2","Gfi1","Ltbr","Klf4","Bcor","Wdr5")
pdf("sel_marker_genes_scatter_all.pdf")
par(mfrow = c(3, 3))
for (ge in sel_mk_genes){
  if (ge %in% rownames(tno_cnt)){
    plot(fate_anno_pDC$bias_cDC1, tno_cnt[ge,], pch = 16,xlab="cDC1 fate bias",ylab="gene expression")
    title(ge)
    plot(fate_anno_pDC$bias_cDC2, tno_cnt[ge,], pch = 16,xlab="cDC2 fate bias",ylab="gene expression")
    title(ge)
    plot(fate_anno_pDC$bias_pDC, tno_cnt[ge,], pch = 16,xlab="pDC fate bias",ylab="gene expression")
    title(ge)
  }
}
dev.off()



pdf("cDC1_genes_bar_all.pdf")
par(mfrow = c(4, 2))
for (ge in ggplot_cDC1_genes){
  plot(fate_anno_pDC$bias_pDC, tno_cnt[ge,], pch = 16,xlab="cDC1 fate bias",ylab="gene expression")
  title(ge)
  #pp = ggplot(data=NULL,aes(x=1:ncol(tno_cnt),y=the_y))+
  #  geom_bar(stat="identity")+
  #  labs(x="cDC1 fate bias",y="gene expression",title=ge)+the_theme
  #print(pp)
  #p_list[[i]] =  pp
  #i = i+1
}
dev.off()

pdf("cDC1_genes_tern_all.pdf")
for (ge in ggplot_cDC1_genes){
  pp = ggtern(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2)) + theme_bw() +
    geom_point(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2,col=tno_cnt[ge,]))+
    scale_color_gradientn(colours = rev(rainbow(5))[2:5])+
    limit_tern(1.05,1.05,1.05)+
    labs(title=ge,x="cDC1",y="pDC",z="cDC2")
  print(pp)
}
dev.off()


pdf("cDC1_genes_tern_all.pdf")
for (ge in ggplot_cDC1_genes){
  pp = ggtern(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2)) + theme_bw() +
    geom_point(data=fate_anno_pDC,alpha=0.85,aes(bias_cDC1,bias_pDC,bias_cDC2,col=tno_cnt[ge,],size=tno_cnt[ge,]))+
    scale_colour_gradient(low="grey",high = "brown4")+
    limit_tern(1.05,1.05,1.05)+
    scale_size_continuous(range=c(2,7))+
    labs(title=ge,x="cDC1",y="pDC",z="cDC2",col="log counts",size="log counts")+
    theme(text = element_text(size=15,face = "bold"))+
    theme_hidegrid()
  print(pp)
}
dev.off()

fate_anno_pDC$pDC_col = ceiling(10*fate_anno_pDC$bias_pDC)/10
fate_anno_pDC$cDC1_col = ceiling(10*fate_anno_pDC$bias_cDC1)/10
fate_anno_pDC$cDC2_col = ceiling(10*fate_anno_pDC$bias_cDC2)/10

fate_anno_pDC$fate_group = paste(fate_anno_pDC$cDC1_col,fate_anno_pDC$cDC2_col,fate_anno_pDC$pDC_col,sep="_")

colo <- rgb(fate_anno_pDC$cDC1_col, fate_anno_pDC$cDC2_col,fate_anno_pDC$pDC_col,alpha=1)



pdf("cDC1_genes_tern_all_tri_color.pdf")
for (ge in ggplot_cDC1_genes){
  pp = ggtern(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2)) + theme_bw() +
    geom_point(data=fate_anno_pDC,alpha=0.7,aes(bias_cDC1,bias_pDC,bias_cDC2,col=fate_group,size=tno_cnt[ge,]))+
    scale_color_manual(guide=FALSE,values = unique(colo),limits = unique(fate_anno_pDC$fate_group))+
    limit_tern(1.05,1.05,1.05)+
    scale_size_continuous(range=c(1,7))+
    labs(title=ge,x="cDC1",y="pDC",z="cDC2",col="log counts",size="log counts")+
    theme(text = element_text(size=15,face = "bold"))+
    theme_hidegrid()
  print(pp)
}
dev.off()





ggplot_cDC2_genes = names(glmout_cDC2)


p_list = list()
i=1

pdf("cDC2_genes_bar_all.pdf")
par(mfrow = c(4, 2))
for (ge in ggplot_cDC2_genes){
  plot(fate_anno_pDC$bias_pDC, tno_cnt[ge,], pch = 16,xlab="cDC2 fate bias",ylab="gene expression")
  title(ge)
  #pp = ggplot(data=NULL,aes(x=1:ncol(tno_cnt),y=the_y))+
  #  geom_bar(stat="identity")+
  #  labs(x="cDC1 fate bias",y="gene expression",title=ge)+the_theme
  #print(pp)
  #p_list[[i]] =  pp
  #i = i+1
}
dev.off()

pdf("cDC2_genes_tern_all.pdf")
for (ge in ggplot_cDC2_genes){
  pp = ggtern(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2)) + theme_bw() +
    geom_point(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2,col=tno_cnt[ge,]))+
    scale_color_gradientn(colours = rev(rainbow(5))[2:5])+
    limit_tern(1.05,1.05,1.05)+
    labs(title=ge,x="cDC1",y="pDC",z="cDC2")
  print(pp)
}
dev.off()


pdf("cDC2_genes_tern_all_tri_color.pdf")
for (ge in ggplot_cDC2_genes){
  pp = ggtern(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2)) + theme_bw() +
    geom_point(data=fate_anno_pDC,alpha=0.7,aes(bias_cDC1,bias_pDC,bias_cDC2,col=fate_group,size=tno_cnt[ge,]))+
    scale_color_manual(guide=FALSE,values = unique(colo),limits = unique(fate_anno_pDC$fate_group))+
    limit_tern(1.05,1.05,1.05)+
    scale_size_continuous(range=c(1,7))+
    labs(title=ge,x="cDC1",y="pDC",z="cDC2",col="log counts",size="log counts")+
    theme(text = element_text(size=15,face = "bold"))+
    theme_hidegrid()
  print(pp)
}
dev.off()






ggplot_pDC_genes = names(glmout_pDC)


p_list = list()
i=1

pdf("pDC_genes_bar_all.pdf")
par(mfrow = c(2, 2))
for (ge in ggplot_pDC_genes){
  #the_y = tno_cnt[ge,order(fate_anno_pDC$bias_pDC)]
  #barplot(unname(the_y),xlab = "cDC2 fate bias")
  
  plot(fate_anno_pDC$bias_pDC, tno_cnt[ge,], pch = 16,xlab="pDC fate bias",ylab="gene expression")
  title(ge)
  #pp = ggplot(data=NULL,aes(x=1:ncol(tno_cnt),y=the_y))+
  #  geom_bar(stat="identity")+
  #  labs(x="cDC1 fate bias",y="gene expression",title=ge)+the_theme
  #print(pp)
  #p_list[[i]] =  pp
  #i = i+1
}
dev.off()

pdf("pDC_genes_tern_all.pdf")
for (ge in ggplot_pDC_genes){
  pp = ggtern(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2)) + theme_bw() +
    geom_point(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2,col=tno_cnt[ge,]))+
    scale_color_gradientn(colours = rev(rainbow(5))[2:5])+
    limit_tern(1.05,1.05,1.05)+
    labs(title=ge,x="cDC1",y="pDC",z="cDC2")
  print(pp)
}
dev.off()



pdf("pDC_genes_tern_all_tri_color.pdf")
for (ge in ggplot_pDC_genes){
  pp = ggtern(data=fate_anno_pDC,aes(bias_cDC1,bias_pDC,bias_cDC2)) + theme_bw() +
    geom_point(data=fate_anno_pDC,alpha=0.7,aes(bias_cDC1,bias_pDC,bias_cDC2,col=fate_group,size=tno_cnt[ge,]))+
    scale_color_manual(guide=FALSE,values = unique(colo),limits = unique(fate_anno_pDC$fate_group))+
    limit_tern(1.05,1.05,1.05)+
    scale_size_continuous(range=c(1,7))+
    labs(title=ge,x="cDC1",y="pDC",z="cDC2",col="log counts",size="log counts")+
    theme(text = element_text(size=15,face = "bold"))+
    theme_hidegrid()
  print(pp)
}
dev.off()



