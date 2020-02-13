test_fc_vec = DE_FC[,c("gene_names","cDC1_logFC")]
test_fc_vec = test_fc_vec[order(test_fc_vec$cDC1_logFC),]

cDC1_gg = selected_genes_3fate$gene_names[grepl("cDC1",selected_genes_3fate$overlap)]

test_fc_vec$cDC1_SISseq = "grey"
test_fc_vec$cDC1_SISseq[test_fc_vec$gene_names %in% cDC1_gg] = "red"
plot(test_fc_vec$cDC1_logFC,col=test_fc_vec$cDC1_SISseq,main="cDC1")



test_fc_vec = DE_FC[,c("gene_names","cDC2_logFC")]
test_fc_vec = test_fc_vec[order(test_fc_vec$cDC2_logFC),]

cDC2_gg = selected_genes_3fate$gene_names[grepl("cDC2",selected_genes_3fate$overlap)]

test_fc_vec$cDC2_SISseq = "grey"
test_fc_vec$cDC2_SISseq[test_fc_vec$gene_names %in% cDC2_gg] = "red"
plot(test_fc_vec$cDC2_logFC,col=test_fc_vec$cDC2_SISseq,main="cDC2")



test_fc_vec = DE_FC[,c("gene_names","cDC1_logFC")]
test_fc_vec = test_fc_vec[order(test_fc_vec$cDC1_logFC),]

cDC1_gg = selected_genes_3fate$gene_names[grepl("cDC1",selected_genes_3fate$overlap)]

test_fc_vec$cDC1_SISseq = "grey"
test_fc_vec$cDC1_SISseq[test_fc_vec$gene_names %in% cDC1_gg] = "red"
plot(test_fc_vec$cDC1_logFC,col=test_fc_vec$cDC1_SISseq,main="cDC1")