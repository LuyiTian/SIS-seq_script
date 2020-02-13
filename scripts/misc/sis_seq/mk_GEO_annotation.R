# make annotation for GEO submission 

load("~/Dropbox/research/sis_seq/clone_split/allcounts.rda")

# including individual sister wells
JT27_all_fate_anno <- read.csv("~/Dropbox/research/sis_seq/fate_anno/JT27/JT27_all_fate_anno.csv")
JS86_all_fate_anno <- read.csv("~/Dropbox/research/sis_seq/fate_anno/JS86/JS86_all_fate_anno.csv")
JS91_all_fate_anno <- read.csv("~/Dropbox/research/sis_seq/fate_anno/JS91/JS91_all_fate_anno.csv")
JS96_all_fate_anno <- read.csv("~/Dropbox/research/sis_seq/fate_anno/JS96/JS96_all_fate_anno.csv")


fc1anno = read.table("~/Dropbox/research/sis_seq/clone_split/r_data/ilsgr_sample_anno Revised.txt", header=TRUE, sep="\t")
fc2anno = read.table("~/Dropbox/research/sis_seq/clone_split/r_data/JS86-DZ11_anno Revised.txt", header=TRUE, sep="\t")
### js91 = read.table("JS91_annotation.txt", header=TRUE, sep="\t")
js91 = read.table("~/Dropbox/research/sis_seq/clone_split/r_data/JS91_sampleanno.txt", header=TRUE, sep="\t") #change from JS91_annotation to JS91_sampleanno
js96 = read.table("~/Dropbox/research/sis_seq/clone_split/r_data/JS96_sampleanno.txt", header=TRUE, sep="\t")



#JT27

#JT27_all_fate_anno$index =  as.factor(rep(1:42,each=2))
#JT27_all_fate_anno$title = paste0("ILSGR0260r_S1_BC",JT27_all_fate_anno$index,".bam")

JT27_all_fate_anno = JT27_all_fate_anno[1:80,] 
# clone index 41 and 42 were not preceed for RNA-seq according to `JS27, Analysis cell numbers-2.xlsx`

JT27_long_anno = data.frame(cDC1_fill_sis1=JT27_all_fate_anno[(1:40)*2-1,"cDC1_fill"],
                            cDC1_fill_sis2=JT27_all_fate_anno[(1:40)*2,"cDC1_fill"],
                            cDC2_fill_sis1=JT27_all_fate_anno[(1:40)*2-1,"cDC2_fill"],
                            cDC2_fill_sis2=JT27_all_fate_anno[(1:40)*2,"cDC2_fill"],
                            pDC_fill_sis1=JT27_all_fate_anno[(1:40)*2-1,"pDC_fill"],
                            pDC_fill_sis2=JT27_all_fate_anno[(1:40)*2,"pDC_fill"],
                            cDC1_sis1=JT27_all_fate_anno[(1:40)*2-1,"cDC1"],
                            cDC1_sis2=JT27_all_fate_anno[(1:40)*2,"cDC1"],
                            cDC2_sis1=JT27_all_fate_anno[(1:40)*2-1,"cDC2"],
                            cDC2_sis2=JT27_all_fate_anno[(1:40)*2,"cDC2"],
                            pDC_sis1=JT27_all_fate_anno[(1:40)*2-1,"pDC"],
                            pDC_sis2=JT27_all_fate_anno[(1:40)*2,"pDC"])

JT27_long_anno$title = paste0("ILSGR0260r_S1_BC",1:40,".bam")

JT27_long_anno$raw_file = paste0("ILSGR0260r_S1_BC",1:40,"_trimmed12_50bp.fastq.gz")

JT27_long_anno$cells_sorted = fc1anno$Number.of.Cells

#JS86
i_max = 39
JS86_long_anno = data.frame(cDC1_fill_sis1=JS86_all_fate_anno[(1:i_max)*2-1,"cDC1_fill"],
                            cDC1_fill_sis2=JS86_all_fate_anno[(1:i_max)*2,"cDC1_fill"],
                            cDC2_fill_sis1=JS86_all_fate_anno[(1:i_max)*2-1,"cDC2_fill"],
                            cDC2_fill_sis2=JS86_all_fate_anno[(1:i_max)*2,"cDC2_fill"],
                            pDC_fill_sis1=JS86_all_fate_anno[(1:i_max)*2-1,"pDC_fill"],
                            pDC_fill_sis2=JS86_all_fate_anno[(1:i_max)*2,"pDC_fill"],
                            cDC1_sis1=JS86_all_fate_anno[(1:i_max)*2-1,"cDC1"],
                            cDC1_sis2=JS86_all_fate_anno[(1:i_max)*2,"cDC1"],
                            cDC2_sis1=JS86_all_fate_anno[(1:i_max)*2-1,"cDC2"],
                            cDC2_sis2=JS86_all_fate_anno[(1:i_max)*2,"cDC2"],
                            pDC_sis1=JS86_all_fate_anno[(1:i_max)*2-1,"pDC"],
                            pDC_sis2=JS86_all_fate_anno[(1:i_max)*2,"pDC"])

JS86_long_anno$title = paste("D9JS27_BC", 1:39, "_R2.bam", sep="")
JS86_long_anno$raw_file = paste0("D9JS27_BC",1:39,"_R2.fastq.gz")

JS86_long_anno$cells_sorted = fc2anno$Number.of.Cells[1:39]
#JS91
JS91_all_fate_anno$index =  as.factor(rep(1:40,each=2))

JS91_long_anno = data.frame(cDC1_fill_sis1=JS91_all_fate_anno[(1:40)*2-1,"cDC1_fill"],
                            cDC1_fill_sis2=JS91_all_fate_anno[(1:40)*2,"cDC1_fill"],
                            cDC2_fill_sis1=JS91_all_fate_anno[(1:40)*2-1,"cDC2_fill"],
                            cDC2_fill_sis2=JS91_all_fate_anno[(1:40)*2,"cDC2_fill"],
                            pDC_fill_sis1=JS91_all_fate_anno[(1:40)*2-1,"pDC_fill"],
                            pDC_fill_sis2=JS91_all_fate_anno[(1:40)*2,"pDC_fill"],
                            cDC1_sis1=JS91_all_fate_anno[(1:40)*2-1,"cDC1"],
                            cDC1_sis2=JS91_all_fate_anno[(1:40)*2,"cDC1"],
                            cDC2_sis1=JS91_all_fate_anno[(1:40)*2-1,"cDC2"],
                            cDC2_sis2=JS91_all_fate_anno[(1:40)*2,"cDC2"],
                            pDC_sis1=JS91_all_fate_anno[(1:40)*2-1,"pDC"],
                            pDC_sis2=JS91_all_fate_anno[(1:40)*2,"pDC"])

JS91_long_anno$title=paste("JS91_", substr(js91[,2],1,6), "_trimmed.bam", sep="")
JS91_long_anno$raw_file = paste0("JS91_",substr(js91[,2],1,6),"_trimmed.fastq.gz")

js91$Cells.Sorted[js91$Cells.Sorted == "blank (neg control)"] = 0
js91$Cells.Sorted = as.numeric(js91$Cells.Sorted)
js91$Cells.Sorted[is.na(js91$Cells.Sorted)] = 0
JS91_long_anno$cells_sorted = js91$Cells.Sorted


#JS96
the_index =  rep(1:39,each=2)
the_index[the_index %in% c(1,2,3,13,21)] = NA
the_index = c(the_index,c(1,1,2,2,NA,NA,3,3,13,13,21,21))
JS96_all_fate_anno$index = the_index
JS96_all_fate_anno = JS96_all_fate_anno[!is.na(JS96_all_fate_anno$index),]

i_max = nrow(JS96_all_fate_anno)/2
JS96_long_anno = data.frame(index = JS96_all_fate_anno$index[(1:i_max)*2-1],
                            cDC1_fill_sis1=JS96_all_fate_anno[(1:i_max)*2-1,"cDC1_fill"],
                            cDC1_fill_sis2=JS96_all_fate_anno[(1:i_max)*2,"cDC1_fill"],
                            cDC2_fill_sis1=JS96_all_fate_anno[(1:i_max)*2-1,"cDC2_fill"],
                            cDC2_fill_sis2=JS96_all_fate_anno[(1:i_max)*2,"cDC2_fill"],
                            pDC_fill_sis1=JS96_all_fate_anno[(1:i_max)*2-1,"pDC_fill"],
                            pDC_fill_sis2=JS96_all_fate_anno[(1:i_max)*2,"pDC_fill"],
                            cDC1_sis1=JS96_all_fate_anno[(1:i_max)*2-1,"cDC1"],
                            cDC1_sis2=JS96_all_fate_anno[(1:i_max)*2,"cDC1"],
                            cDC2_sis1=JS96_all_fate_anno[(1:i_max)*2-1,"cDC2"],
                            cDC2_sis2=JS96_all_fate_anno[(1:i_max)*2,"cDC2"],
                            pDC_sis1=JS96_all_fate_anno[(1:i_max)*2-1,"pDC"],
                            pDC_sis2=JS96_all_fate_anno[(1:i_max)*2,"pDC"])

JS96_long_anno = JS96_long_anno[order(JS96_long_anno$index),]

JS96_long_anno = JS96_long_anno[,!(colnames(JS96_long_anno) == "index")]

JS96_long_anno$title = paste("JS96_", substr(js96[1:39,2],1,6), "_trimmed.bam", sep="") # index 40 is blank control, have no data
JS96_long_anno$raw_file = paste0("JS96_",substr(js96[1:39,2],1,6),"_trimmed.fastq.gz")

js96$Cells.Sorted[js96$Cells.Sorted == "None"] = 0
js96$Cells.Sorted = as.numeric(js96$Cells.Sorted)
js96$Cells.Sorted[is.na(js96$Cells.Sorted)] = 0
JS96_long_anno$cells_sorted = js96$Cells.Sorted[1:39]

JT27_long_anno$batch = "JT27"
JS86_long_anno$batch = "JS86"
JS91_long_anno$batch = "JS91"
JS96_long_anno$batch = "JS96"

conbined_long_anno = rbind(JT27_long_anno,JS86_long_anno,JS91_long_anno,JS96_long_anno)

conbined_long_anno$title = unlist(lapply(strsplit(conbined_long_anno$title,"[.]"),function(x){x[1]}))
conbined_long_anno$Sample_name = conbined_long_anno$title
conbined_long_anno$source_name = "hematopoietic stem and progenitor cells"
conbined_long_anno$organism = "Mus musculus"
conbined_long_anno$molecule = "polyA RNA"
conbined_long_anno$description = "progenitor cells at day 2.5"
conbined_long_anno$description[conbined_long_anno$cells_sorted==0] = "no cell sorted"
conbined_long_anno$processed_data_file = "gene_counts.csv"
conbined_long_anno = conbined_long_anno[,c("Sample_name","title","source_name","organism",
                                           "cDC1_fill_sis1","cDC1_fill_sis2","cDC2_fill_sis1","cDC2_fill_sis2","pDC_fill_sis1","pDC_fill_sis2",
                                           "cDC1_sis1","cDC1_sis2","cDC2_sis1","cDC2_sis2","pDC_sis1","pDC_sis2", "cells_sorted",
                                           "molecule","description","processed_data_file","raw_file")]

write.csv(conbined_long_anno,file="meta_anno.csv",row.names = F)

fc_sel = fc
colnames(fc_sel) = unlist(lapply(strsplit(colnames(fc_sel),"[.]"),function(x){x[1]}))

fc_sel = fc_sel[,colnames(fc_sel) %in% conbined_long_anno$Sample_name]
fc_sel = fc_sel[,conbined_long_anno$Sample_name]
write.csv(fc_sel,file="gene_counts.csv")

library(ggplot2)
library(ggpubr)

p1 = ggplot(data=conbined_long_anno,aes(x=cDC1_sis1+1,y=cDC1_sis2+1,col=batch))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="sister well a",y="sister well b",title="normalized cDC1 cell counts")+
  theme(text = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2 = ggplot(data=conbined_long_anno,aes(x=cDC2_sis1+1,y=cDC2_sis2+1,col=batch))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="sister well a",y="sister well b",title="normalized cDC2 cell counts")+
  theme(text = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p3 = ggplot(data=conbined_long_anno,aes(x=pDC_sis1+1,y=pDC_sis2+1,col=batch))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  labs(x="sister well a",y="sister well b",title="normalized pDC cell counts")+
  theme(text = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#pdf("supp_fig1.pdf",width = 4,height = 10)
#ggarrange(p1, p2, p3, 
#          labels = c("A", "B", "C"),
#         ncol = 1, nrow = 3)
#dev.off()



# conbined_long_anno_tmp = conbined_long_anno
# 
# conbined_long_anno_tmp$sis1_total = conbined_long_anno_tmp$cDC1_sis1+conbined_long_anno_tmp$cDC2_sis1+conbined_long_anno_tmp$pDC_sis1
# conbined_long_anno_tmp$sis2_total = conbined_long_anno_tmp$cDC1_sis2+conbined_long_anno_tmp$cDC2_sis2+conbined_long_anno_tmp$pDC_sis2
# 
# conbined_long_anno_tmp = conbined_long_anno_tmp[conbined_long_anno_tmp$sis1_total+conbined_long_anno_tmp$sis2_total>1,]
# 
# pdf("supp_fig1.pdf",width = 5.5,height = 5)
# ggplot(data=conbined_long_anno_tmp,aes(x=sis1_total+1,y=sis2_total+1,col=batch))+
#   geom_point()+
#   scale_x_log10()+
#   scale_y_log10()+
#   labs(x="sister well a",y="sister well b",title="normalized total cell counts")+
#   theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# dev.off()



checklist.md5 <- read.table("~/Dropbox/research/sis_seq/checklist.md5.tsv", quote="\"", comment.char="")

checklist.md5[!(checklist.md5$V2 %in% conbined_long_anno$raw_file),]
rownames(checklist.md5) = checklist.md5$V2
checklist.md5 = checklist.md5[(checklist.md5$V2 %in% conbined_long_anno$raw_file),]
checklist.md5 = checklist.md5[conbined_long_anno$raw_file,]

fq_info = data.frame(file_name=conbined_long_anno$raw_file,file_type=rep("fastq",nrow(conbined_long_anno)),checksum=checklist.md5$V1,inst="NextSeq 500",stringsAsFactors = F)

fq_info$inst[grepl("ILSGR0260r", fq_info$file_name)] = "MiSeq"

write.csv(fq_info,file="fq_anno.csv",row.names = F)
