# 
# cDC1_genes <- read.table("~/Dropbox/research/sis_seq/cDC1_genes.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)
# cDC2_genes <- read.table("~/Dropbox/research/sis_seq/cDC2_genes.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)
# pDC_genes <- read.table("~/Dropbox/research/sis_seq/pDC_genes.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)
# 
# SIS_seq_genes = unique(c(cDC1_genes$V1,
#                       cDC2_genes$V1,
#                       pDC_genes$V1))
# write.table(SIS_seq_genes,file="SIS_seq_genes.csv",row.names = F, quote = F,col.names = F)
# 
cDC1_Ig_genes <- read.table("~/Dropbox/research/sis_seq/cDC1_Ig_genes.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)
cDC2_Ig_genes <- read.table("~/Dropbox/research/sis_seq/cDC2_Ig_genes.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)
pDC_Ig_genes <- read.table("~/Dropbox/research/sis_seq/pDC_Ig_genes.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)
# 
CDP_genes_Ig <- read.table("~/Dropbox/research/sis_seq/CDP_genes_Ig.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)

# 
Nat_Imm_genes <- read.table("~/Dropbox/research/sis_seq/Nat_Imm_genes.csv", quote="\"", comment.char="", stringsAsFactors=FALSE)

Ig_mature_genes = unique(c(cDC1_Ig_genes$V1,
                  cDC2_Ig_genes$V1,
                  pDC_Ig_genes$V1))
write.table(Ig_mature_genes,file="Ig_mature_genes.csv",row.names = F, quote = F,col.names = F)


selected_genes_3fate <- read.csv("~/Dropbox/research/sis_seq/selected_genes_3fate.csv", stringsAsFactors=FALSE)

SIS_seq_genes= unique(selected_genes_3fate$gene_names[(grepl("cDC1_SISeq", selected_genes_3fate$overlap) | grepl("cDC2_SISeq", selected_genes_3fate$overlap)) | grepl("pDC_SISeq", selected_genes_3fate$overlap)])

write.table(SIS_seq_genes,file="SIS_seq_genes.csv",row.names = F, quote = F,col.names = F)


library(gplots)
venn(list(SIS_seq_genes=SIS_seq_genes,
            Blood=coexpressed_Blood,
            Nat=coexpressed_Nat),intersection=TRUE)

write.csv(SIS_seq_genes,file="SIS_seq_genes.csv",row.names = FALSE,quote=FALSE)
write.csv(coexpressed_Blood,file="coexpressed_Blood.csv",row.names = FALSE,quote=FALSE)
write.csv(coexpressed_Nat,file="coexpressed_Nat.csv",row.names = FALSE,quote=FALSE)
