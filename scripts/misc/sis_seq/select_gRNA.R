library(biomaRt)

ensembl = useEnsembl(biomart="ensembl",dataset="mmusculus_gene_ensembl")

sele_genes <- getBM(attributes=c('ensembl_gene_id',
                                 'mgi_id','external_gene_name','ensembl_exon_id'), filters =
                      'external_gene_name', values =unique(c(res_df$gene_names,interest_genes)), mart = ensembl)

Mouse_Sanger_Informatics <- read_excel("Dropbox/research/sis_seq/Mouse Sanger Informatics.xlsx")

Mouse_sel = Mouse_Sanger_Informatics[Mouse_Sanger_Informatics$`ensembl gene_id` %in% sele_genes$mgi_id,]

length(unique(Mouse_sel$`ensembl gene_id`))
length(unique(sele_genes$mgi_id))
length(unique(c(res_df$gene_names,interest_genes)))

unique(sele_genes$external_gene_name[!(sele_genes$mgi_id %in% Mouse_Sanger_Informatics$`ensembl gene_id`)])
# [1] "Snord22"       "2610307P16Rik" "Akap5"         "Ggn"           "Zeb2"          "Gpank1"        "Zscan29"       "Rag2"         
# [9] "Ccr9"          "Spsb1"         "Ccr2"          "Crem"          "Pde4dip"       "Tsc22d1"       "Zfp263"        "Cebpd"        
# [17] "Zfp932"        "Zfp579"        "Cav1"          "Pstpip2"   


Mouse_sel_1g = Mouse_sel[!duplicated(Mouse_sel$`ensembl gene_id`),]

write.csv(Mouse_sel_1g, file="mouse_select_189gRNA.csv",row.names = F)

#####
sele_genes <- getBM(attributes=c('ensembl_gene_id',
                                 'mgi_id','external_gene_name','ensembl_exon_id'), filters =
                      'external_gene_name', values =unique(c(res_df0$gene_names,interest_genes)), mart = ensembl)

Mouse_Sanger_Informatics <- read_excel("Dropbox/research/sis_seq/Mouse Sanger Informatics.xlsx")

Mouse_sel = Mouse_Sanger_Informatics[Mouse_Sanger_Informatics$`ensembl gene_id` %in% sele_genes$mgi_id,]


Mouse_sel_1g = Mouse_sel[!duplicated(Mouse_sel$`ensembl gene_id`),]

write.csv(Mouse_sel_1g, file="mouse_select_445gRNA.csv",row.names = F)

write.csv(Mouse_sel, file="mouse_select_445X2_gRNA.csv",row.names = F)
