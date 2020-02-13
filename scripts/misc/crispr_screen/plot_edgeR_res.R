#plot edgeR DE
library(ggplot2)
library(edgeR)
library(ggrepel)
library(ggpubr)
#setwd("~/Dropbox/research/sis_seq/validation/CRISPR_screen")
#sel_mk_genes = c("Cbfa2t3","Mycl2","Id2","Batf3","Irf8","Runx1","Zbtb46","Runx2","Tcf4","Zeb2","Gfi1","Ltbr","Klf4","Bcor","Wdr5")

#known_marker = c("Runx2","Runx1","Gfi1","Irf8","Notch2","Relb","Tcf4","Zeb2","Brd3","Zbtb46","Pqlc3","Plac8","Jade1","Spib")

sel_mk_genes = c("Zc3h12a","Ccl2","Cbfa2t3","Mycl","Id2","Batf3","Irf8","Runx1","Zbtb46","Runx2","Tcf4","Zeb2","Gfi1","Ltbr","Klf4","Bcor","Wdr5")

sis_seq_known = known_list[known_list %in% sis_seq_all]
internal_control = known_list[!(known_list %in% sis_seq_all)]

plot_df = tp1$table
plot_df = plot_df[plot_df$logCPM>5,]
plot_df$nlp = -log10(plot_df$FDR)
plot_df$nlp[plot_df$nlp>80] = 80
plot_df$known_control = "SIS-seq (novel)"
plot_df$known_control[plot_df$gene_names %in% sis_seq_known] = "SIS-seq (known)"
plot_df$known_control[plot_df$gene_names %in% internal_control] = "Internal control"
plot_df$sel_genes = plot_df$gene_names
plot_df$sel_genes[!(plot_df$sel_genes %in% sel_mk_genes) & plot_df$known_control=="SIS-seq (novel)"] = ""
plot_df$sel_genes[abs(plot_df$logFC)<1.5] = ""
plot_df$sel_genes[-log10(plot_df$FDR)<20] = ""
#plot_df$sel_genes[-log10(plot_df$FDR)>40] = plot_df$gene_names[-log10(plot_df$FDR)>40]

p1 = ggplot(data=plot_df,aes(x=logFC,y=nlp,col=known_control,label=sel_genes))+
  geom_point(alpha=0.8,size=1,show.legend = TRUE)+
  geom_text_repel(show.legend = FALSE)+
  labs(y="-log10(FDR)",col="")+
  scale_color_manual(values = c("#1B9E77","#DA9389","#000000"))+
  ylim(0,80)+
  ggtitle("cDC1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1

plot_df = tp2$table
plot_df = plot_df[plot_df$logCPM>5,]
plot_df$nlp = -log10(plot_df$FDR)
plot_df$nlp[plot_df$nlp>80] = 80
plot_df$known_control = "SIS-seq (novel)"
plot_df$known_control[plot_df$gene_names %in% sis_seq_known] = "SIS-seq (known)"
plot_df$known_control[plot_df$gene_names %in% internal_control] = "Internal control"
plot_df$sel_genes = plot_df$gene_names
plot_df$sel_genes[!(plot_df$sel_genes %in% sel_mk_genes) & plot_df$known_control=="SIS-seq (novel)"] = ""
plot_df$sel_genes[abs(plot_df$logFC)<1.5] = ""
plot_df$sel_genes[-log10(plot_df$FDR)<20] = ""
#plot_df$sel_genes[-log10(plot_df$FDR)>50] = plot_df$gene_names[-log10(plot_df$FDR)>50]

p2 = ggplot(data=plot_df,aes(x=logFC,y=nlp,col=known_control,label=sel_genes))+
  geom_point(alpha=0.8,size=1,show.legend = TRUE)+
  geom_text_repel(show.legend = FALSE)+
  labs(y="-log10(FDR)",col="")+
  scale_color_manual(values = c("#1B9E77","#DA9389","#000000"))+
  ylim(0,80)+
  ggtitle("cDC2")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2

plot_df = tp3$table
plot_df = plot_df[plot_df$logCPM>5,]
plot_df$nlp = -log10(plot_df$FDR)
plot_df$nlp[plot_df$nlp>80] = 80
plot_df$known_control = "SIS-seq (novel)"
plot_df$known_control[plot_df$gene_names %in% sis_seq_known] = "SIS-seq (known)"
plot_df$known_control[plot_df$gene_names %in% internal_control] = "Internal control"
plot_df$sel_genes = plot_df$gene_names
plot_df$sel_genes[!(plot_df$sel_genes %in% sel_mk_genes) & plot_df$known_control=="SIS-seq (novel)"] = ""
plot_df$sel_genes[abs(plot_df$logFC)<1.5] = ""
plot_df$sel_genes[-log10(plot_df$FDR)<20] = ""
#plot_df$sel_genes[-log10(plot_df$FDR)>60] = plot_df$gene_names[-log10(plot_df$FDR)>60]

p3 = ggplot(data=plot_df,aes(x=logFC,y=nlp,col=known_control,label=sel_genes))+
  geom_point(alpha=0.8,size=1,show.legend = TRUE)+
  geom_text_repel(show.legend = FALSE)+
  labs(y="-log10(FDR)",col="")+
  scale_color_manual(values = c("#1B9E77","#DA9389","#000000"))+
  ylim(0,80)+
  ggtitle("pDC")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3

ggarrange(p1,p2,p3,ncol=3,nrow=1,common.legend = TRUE,align = "v")
ggsave(filename="/Users/tian.l/Dropbox/research/sis_seq/validation/CRISPR_screen/figs/edgeR_volcano_plots.pdf",width = 8.5,height = 4.5)





