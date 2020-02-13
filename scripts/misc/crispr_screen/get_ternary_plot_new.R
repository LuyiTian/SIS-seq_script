library(ggtern)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(ggpubr)
library(plyr)
library(tidyr)
library(RColorBrewer)
hist(log2(rowSums(mtx[,c(1,2,5,6,9,10)])))
hist(log2(rowSums(mtx[,c(1,2,5,6,9,10)+2])))
# choose 2^8 as cutoff

#kp1 = rowSums(mtx[,c(1,2,5,6,9,10)]) > 2^9
#kp2 = rowSums(mtx[,c(1,2,5,6,9,10)+1]) > 2^9

g_na = dge$genes$gene_names#[kp1 & kp2]

#prop_diff = as.data.frame(mtx[kp1& kp2,-ncol(mtx)])
prop_diff = as.data.frame(mtx[,-ncol(mtx)])

prop_diff$cDC1_WT = prop_diff$WT.cDC1.A+prop_diff$WT.cDC1.B
prop_diff$cDC2_WT = prop_diff$WT.cDC2.A+prop_diff$WT.cDC2.B
prop_diff$pDC_WT = prop_diff$WT.pDC.A+prop_diff$WT.pDC.B
prop_diff$WT_sum = rowSums(prop_diff[,c("cDC1_WT","cDC2_WT","pDC_WT")])
prop_diff$cDC1_WT = prop_diff$cDC1_WT/prop_diff$WT_sum
prop_diff$cDC2_WT = prop_diff$cDC2_WT/prop_diff$WT_sum
prop_diff$pDC_WT = prop_diff$pDC_WT/prop_diff$WT_sum

prop_diff$cDC1_Cas9 = prop_diff$Cas9.cDC1.A+prop_diff$Cas9.cDC1.B
prop_diff$cDC2_Cas9 = prop_diff$Cas9.cDC2.A+prop_diff$Cas9.cDC2.B
prop_diff$pDC_Cas9 = prop_diff$Cas9.pDC.A+prop_diff$Cas9.pDC.B
prop_diff$Cas9_sum = rowSums(prop_diff[,c("cDC1_Cas9","cDC2_Cas9","pDC_Cas9")])
prop_diff$cDC1_Cas9 = prop_diff$cDC1_Cas9/prop_diff$Cas9_sum
prop_diff$cDC2_Cas9 = prop_diff$cDC2_Cas9/prop_diff$Cas9_sum
prop_diff$pDC_Cas9 = prop_diff$pDC_Cas9/prop_diff$Cas9_sum

prop_diff$cDC1_diff = prop_diff$cDC1_Cas9-prop_diff$cDC1_WT
prop_diff$cDC2_diff = prop_diff$cDC2_Cas9-prop_diff$cDC2_WT
prop_diff$pDC_diff = prop_diff$pDC_Cas9-prop_diff$pDC_WT

#######

known_marker = c("Runx2","Runx1","Gfi1","Irf8","Notch2","Relb","Tcf4","Zeb2","Brd3","Zbtb46","Pqlc3","Plac8","Jade1","Spib")
sel_mk_genes = c("Zc3h12a","Ccl2","Cbfa2t3","Mycl","Id2","Batf3","Irf8","Runx1","Zbtb46","Runx2","Tcf4","Zeb2","Gfi1","Ltbr","Klf4","Bcor","Wdr5")

prop_diff_ge = prop_diff[prop_diff$Cas9_sum>150,c("cDC1_Cas9", "cDC2_Cas9", "pDC_Cas9")]
prop_diff_ge$gene_name = g_na[prop_diff$Cas9_sum>150]
prop_diff_ge$gene_name[!(prop_diff_ge$gene_name %in% c(known_marker,sel_mk_genes))] = ""
prop_diff_ge$gene_name[!(apply(prop_diff_ge[,c("cDC1_Cas9", "cDC2_Cas9", "pDC_Cas9")],1,min)<0.2 | apply(prop_diff_ge[,c("cDC1_Cas9", "cDC2_Cas9", "pDC_Cas9")],1,max)>0.65)] = ""
prop_diff_ge$sel_ge = "not"
prop_diff_ge$sel_ge[!(prop_diff_ge$gene_name == "")] = "sel"

prop_diff_kn = prop_diff_ge[!(prop_diff_ge$gene_name == ""),]
#prop_diff_kn = prop_diff_kn %>% 
#  group_by(gene_name) %>% 
#  summarise_at(c("cDC1_Cas9", "cDC2_Cas9", "pDC_Cas9"), mean)
prop_diff_kn$known_control = "No"
prop_diff_kn$known_control[prop_diff_kn$gene_name %in% known_marker] = "Yes"

prop_diff_wt = prop_diff[prop_diff$WT_sum>150,c("cDC1_WT","cDC2_WT","pDC_WT")]

cols <- c("not" = "grey", "sel" = "black", "No" = brewer.pal(3, "Dark2")[2], "Yes" = brewer.pal(3, "Dark2")[1])

ggtern()+
  geom_point(data=prop_diff_ge,aes(x=cDC1_Cas9,y=pDC_Cas9,z=cDC2_Cas9,col=sel_ge),alpha=0.8,size=1,show.legend = FALSE)+
  geom_text(data=prop_diff_kn,aes(x=cDC1_Cas9, y=pDC_Cas9, z=cDC2_Cas9,col=known_control,label=gene_name),show.legend = TRUE)+
  stat_density_tern(data=prop_diff_wt,aes(x=cDC1_WT,y=pDC_WT,z=cDC2_WT),alpha=0.7)+
  scale_color_manual(values = cols)+
  scale_fill_brewer(palette = "Dark2",direction = -1)+
  theme_bw()+
  theme_hidegrid()
ggsave(filename="/Users/tian.l/Dropbox/research/sis_seq/validation/CRISPR_screen/figs/ternary_plots_gRNA_new.pdf",width = 5,height = 5)

#######

sel_genes = c("Zfp455","Ifi205","Zscan29","Tyrobp","Klf12","Ankrd24","Notch2","Nrxn1",
              "Zc3h12a","Adam19","Naaa","Trpa1","Celsr3","Mbd5","Ubr3","Drd1","Impdh2","Zfp652","Zfp831")

sel_snd = c("Havcr2","F2rl1","Il12rb1","Csf1r","Eif4ebp3","Tmem108","Wnt2b","Plxna2","Siglecg","Tespa1","Cfc1","Llgl1","Ccdc150")

data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#####
sel_gRNA = c("737","740","121","122","124","470","471","472","351","352","2237",
             "2238","2239","2240","345","346","347","348","113","114","115","805","806","962","964","433","435","436",
             "1377","1378","1379","1380","1713","1714","1716","430","431","2321","2322","893","894","896","1914","1916")
cnt_diff = as.data.frame(mtx[,-ncol(mtx)])

cnt_diff$cDC1_Cas9 = cnt_diff$Cas9.cDC1.A+cnt_diff$Cas9.cDC1.B
cnt_diff$cDC2_Cas9 = cnt_diff$Cas9.cDC2.A+cnt_diff$Cas9.cDC2.B
cnt_diff$pDC_Cas9 = cnt_diff$Cas9.pDC.A+cnt_diff$Cas9.pDC.B
cnt_diff$cDC1_WT = cnt_diff$WT.cDC1.A+cnt_diff$WT.cDC1.B
cnt_diff$cDC2_WT = cnt_diff$WT.cDC2.A+cnt_diff$WT.cDC2.B
cnt_diff$pDC_WT = cnt_diff$WT.pDC.A+cnt_diff$WT.pDC.B

ix = 1
p_list = list()

for (a_g in known_marker){
  cnt_diff_ag = cnt_diff[g_na==a_g,]
  cnt_diff_ag = cnt_diff_ag[,c("cDC1_Cas9","cDC2_Cas9","pDC_Cas9","cDC1_WT","cDC2_WT","pDC_WT")]
  cnt_diff_ag = cnt_diff_ag[rowSums(cnt_diff_ag)>800,]
  cnt_diff_ag = cnt_diff_ag[rownames(cnt_diff_ag) %in% sel_gRNA,]
  if (nrow(cnt_diff_ag)>0){
    nr = nrow(cnt_diff_ag)
    cnt_diff_ag$cDC1_Cas9 = cnt_diff_ag$cDC1_Cas9/cnt_diff_ag$cDC1_WT
    cnt_diff_ag$cDC1_WT = cnt_diff_ag$cDC1_WT/(cnt_diff_ag$cDC1_WT)
    cnt_diff_ag$cDC2_Cas9 = cnt_diff_ag$cDC2_Cas9/cnt_diff_ag$cDC2_WT
    cnt_diff_ag$cDC2_WT = cnt_diff_ag$cDC2_WT/(cnt_diff_ag$cDC2_WT)
    cnt_diff_ag$pDC_Cas9 = cnt_diff_ag$pDC_Cas9/cnt_diff_ag$pDC_WT
    cnt_diff_ag$pDC_WT = cnt_diff_ag$pDC_WT/(cnt_diff_ag$pDC_WT)
    
    
    cnt_diff_long_r <- gather(cnt_diff_ag, condition, gRNA_ratio, cDC1_Cas9:pDC_WT, factor_key=TRUE)
    cnt_diff_long <- data_summary(cnt_diff_long_r, varname="gRNA_ratio", 
                                  groupnames=c("condition"))
    cnt_diff_long = separate(data = cnt_diff_long, col = condition, into = c("cell_type", "cond"),sep="_")
    cnt_diff_long_r = separate(data = cnt_diff_long_r, col = condition, into = c("cell_type", "cond"),sep="_")
    #cnt_diff_long$cond <- ordered(cnt_diff_long$cond, levels = c("cDC1",  "cDC2","pDC"))
    
    p_list[[ix]] = ggplot(data=cnt_diff_long,aes(x=cell_type,y=gRNA_ratio,fill=cond))+
      geom_bar(alpha=0.7,stat="identity", color="black", 
               position=position_dodge()) +
      geom_point(data=cnt_diff_long_r,aes(x=cell_type,y=gRNA_ratio,col=cond), position = position_jitterdodge(jitter.width=0.1),size=0.7,show.legend = F)+
      geom_errorbar(aes(ymin=gRNA_ratio-sd, ymax=gRNA_ratio+sd), width=.2,
                    position=position_dodge(.9)) +
      labs(fill="conditions:")+
      ggtitle(paste0(a_g," (n=",nr,")"))+
      #scale_fill_brewer(palette = "Set1")+
      scale_fill_manual(values = c("#D95F02","#613CB3"))+
      scale_color_manual(values = c("black","black"))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),axis.title.y=element_blank())
    ix = ix+1
  }
}
#ggarrange(plotlist = p_list)

ix = 1
p_list_other = list()
#other_genes = unique(c(sel_genes, sel_snd, sel_mk_genes))
#other_genes = other_genes[!(other_genes %in% known_marker)]
other_genes = c("Id2","Batf3","Cbfa2t3","Klf4","Zfp831","Wnt2b","Bcor","Zc3h12a")
for (a_g in other_genes){
  cnt_diff_ag = cnt_diff[g_na==a_g,]
  cnt_diff_ag = cnt_diff_ag[,c("cDC1_Cas9","cDC2_Cas9","pDC_Cas9","cDC1_WT","cDC2_WT","pDC_WT")]
  cnt_diff_ag = cnt_diff_ag[rowSums(cnt_diff_ag)>800,]
  cnt_diff_ag = cnt_diff_ag[rownames(cnt_diff_ag) %in% sel_gRNA,]
  if (nrow(cnt_diff_ag)>0){
    nr = nrow(cnt_diff_ag)
    cnt_diff_ag$cDC1_Cas9 = cnt_diff_ag$cDC1_Cas9/cnt_diff_ag$cDC1_WT
    cnt_diff_ag$cDC1_WT = cnt_diff_ag$cDC1_WT/(cnt_diff_ag$cDC1_WT)
    cnt_diff_ag$cDC2_Cas9 = cnt_diff_ag$cDC2_Cas9/cnt_diff_ag$cDC2_WT
    cnt_diff_ag$cDC2_WT = cnt_diff_ag$cDC2_WT/(cnt_diff_ag$cDC2_WT)
    cnt_diff_ag$pDC_Cas9 = cnt_diff_ag$pDC_Cas9/cnt_diff_ag$pDC_WT
    cnt_diff_ag$pDC_WT = cnt_diff_ag$pDC_WT/(cnt_diff_ag$pDC_WT)
    
    
    cnt_diff_long_r <- gather(cnt_diff_ag, condition, gRNA_ratio, cDC1_Cas9:pDC_WT, factor_key=TRUE)
    cnt_diff_long <- data_summary(cnt_diff_long_r, varname="gRNA_ratio", 
                                  groupnames=c("condition"))
    cnt_diff_long = separate(data = cnt_diff_long, col = condition, into = c("cell_type", "cond"),sep="_")
    cnt_diff_long_r = separate(data = cnt_diff_long_r, col = condition, into = c("cell_type", "cond"),sep="_")
    #cnt_diff_long$cond <- ordered(cnt_diff_long$cond, levels = c("cDC1",  "cDC2","pDC"))
    
    p_list_other[[ix]] = ggplot(data=cnt_diff_long,aes(x=cell_type,y=gRNA_ratio,fill=cond))+
      geom_bar(alpha=0.7,stat="identity", color="black", 
               position=position_dodge()) +
      geom_point(data=cnt_diff_long_r,aes(x=cell_type,y=gRNA_ratio,col=cond), position = position_jitterdodge(jitter.width=0.1),size=0.7,show.legend = F)+
      geom_errorbar(aes(ymin=gRNA_ratio-sd, ymax=gRNA_ratio+sd), width=.2,
                    position=position_dodge(.9)) +
      labs(fill="conditions:")+
      ggtitle(paste0(a_g," (n=",nr,")"))+
      #scale_fill_brewer(palette = "Set1")+
      scale_fill_manual(values = c("#D95F02","#613CB3"))+
      scale_color_manual(values = c("black","black"))+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), panel.border = element_blank(),axis.line = element_line(colour = "black"),
            axis.title.x=element_blank(),axis.title.y=element_blank())
    ix = ix+1
  }
}
ggarrange(plotlist = append(p_list,p_list_other),ncol=4,nrow=3,align="v",common.legend=TRUE)

ggsave(filename="/Users/tian.l/Dropbox/research/sis_seq/validation/CRISPR_screen/figs/individual_plots_gRNA1.pdf",width = 7,height = 6)
#####

###### for Sara
pdf("for_sara_individual_plots.pdf")
for (a_g in c(sel_second_batch,known_marker,sel_mk_genes,"Batf3","Irf8","Id2","Bcor")){
  cnt_diff_ag = cnt_diff[g_na==a_g,]
  cnt_diff_ag = cnt_diff_ag[,c("cDC1_Cas9","cDC2_Cas9","pDC_Cas9","cDC1_WT","cDC2_WT","pDC_WT")]
  cnt_diff_ag = cnt_diff_ag[rowSums(cnt_diff_ag)>800,]
  #cnt_diff_ag = cnt_diff_ag[rownames(cnt_diff_ag) %in% sel_gRNA,]
  if (nrow(cnt_diff_ag)>0){
    nr = nrow(cnt_diff_ag)
    cnt_diff_ag$cDC1_Cas9 = cnt_diff_ag$cDC1_Cas9/cnt_diff_ag$cDC1_WT
    cnt_diff_ag$cDC1_WT = cnt_diff_ag$cDC1_WT/(cnt_diff_ag$cDC1_WT)
    cnt_diff_ag$cDC2_Cas9 = cnt_diff_ag$cDC2_Cas9/cnt_diff_ag$cDC2_WT
    cnt_diff_ag$cDC2_WT = cnt_diff_ag$cDC2_WT/(cnt_diff_ag$cDC2_WT)
    cnt_diff_ag$pDC_Cas9 = cnt_diff_ag$pDC_Cas9/cnt_diff_ag$pDC_WT
    cnt_diff_ag$pDC_WT = cnt_diff_ag$pDC_WT/(cnt_diff_ag$pDC_WT)
    
    
    cnt_diff_long <- gather(cnt_diff_ag, condition, gRNA_ratio, cDC1_Cas9:pDC_WT, factor_key=TRUE)
    cnt_diff_long <- data_summary(cnt_diff_long, varname="gRNA_ratio", 
                                  groupnames=c("condition"))
    cnt_diff_long = separate(data = cnt_diff_long, col = condition, into = c("cell_type", "cond"),sep="_")
    #cnt_diff_long$cond <- ordered(cnt_diff_long$cond, levels = c("cDC1",  "cDC2","pDC"))
    
    p = ggplot(data=cnt_diff_long,aes(x=cell_type,y=gRNA_ratio,fill=cond))+
      geom_bar(alpha=0.7,stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=gRNA_ratio-sd, ymax=gRNA_ratio+sd), width=.2,
                    position=position_dodge(.9)) +
      labs(fill="conditions:")+
      ggtitle(paste0(a_g," (n=",nr,")"))+
      scale_fill_brewer(palette = "Set1")+
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())
    print(p)
    ix = ix+1
  }
}
dev.off()
######



pairs(prop_diff[,c("cDC1_diff","cDC2_diff","pDC_diff")])

null_mtx = prop_diff[,c("cDC1_WT", "cDC2_WT", "pDC_WT")]

null_mtx_r = null_mtx[sample(rownames(null_mtx)),]

null_mtx = null_mtx_r-null_mtx

WT_sig_val = list()
WT_sig_val$cDC1_lo = unname(quantile(null_mtx$cDC1_WT,0.05))
WT_sig_val$cDC1_hi = unname(quantile(null_mtx$cDC1_WT,0.95))
WT_sig_val$cDC2_lo = unname(quantile(null_mtx$cDC2_WT,0.05))
WT_sig_val$cDC2_hi = unname(quantile(null_mtx$cDC2_WT,0.95))
WT_sig_val$pDC_lo =  unname(quantile(null_mtx$pDC_WT,0.05))
WT_sig_val$pDC_hi =  unname(quantile(null_mtx$pDC_WT,0.95))

get_cosine_dist = function(prop_mtx){
  X = prop_mtx
  cos.sim <- function(ix) 
  {
    
    A = X[ix[1],]
    B = X[ix[2],]
    if (sqrt(sum(A^2))<0.1 | sqrt(sum(B^2))<0.1){return(0)}
    return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
  }   
  n <- nrow(X) 
  cmb <- expand.grid(i=1:n, j=1:n) 
  C <- matrix(apply(cmb,1,cos.sim),n,n)
  C[upper.tri(C)]
}

all_dist = c()
gRNA_fate_df = NULL
for (i in unique(g_na)){
  tmp_v = rep(FALSE,6) 
  prop_mtx=prop_diff[g_na == i,c("cDC1_diff","cDC2_diff","pDC_diff")]
  di_tmp = get_cosine_dist(prop_mtx)
  all_dist = c(all_dist,max(di_tmp))
  if(max(di_tmp)<0.8){
    next
  }
  cDC1_lo = unname(table(prop_mtx[,1] < WT_sig_val$cDC1_lo)["TRUE"])
  cDC1_hi = unname(table(prop_mtx[,1] > WT_sig_val$cDC1_hi)["TRUE"])
  
  if (!is.na(cDC1_lo)){
    if (cDC1_lo>1){tmp_v[1] = TRUE}
  }else if(!is.na(cDC1_hi)){
    if (cDC1_hi>1){tmp_v[2] = TRUE}
  }
  
  cDC2_lo = unname(table(prop_mtx[,2] < WT_sig_val$cDC2_lo)["TRUE"])
  cDC2_hi = unname(table(prop_mtx[,2] > WT_sig_val$cDC2_hi)["TRUE"])
  if (!is.na(cDC2_lo)){
    if (cDC2_lo>1){tmp_v[3] = TRUE}
  }else if(!is.na(cDC2_hi)){
    if (cDC2_hi>1){tmp_v[4] = TRUE}
  }
  
  pDC_lo = unname(table(prop_mtx[,3] < WT_sig_val$pDC_lo)["TRUE"])
  pDC_hi = unname(table(prop_mtx[,3] > WT_sig_val$pDC_hi)["TRUE"])
  if (!is.na(pDC_lo)){
    if (pDC_lo>1){tmp_v[5] = TRUE}
  }else if(!is.na(pDC_hi)){
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
gRNA_fate_df$keep = apply(gRNA_fate_df,1,function(x){table(x[2:7])["FALSE"]<6})
gRNA_fate_df = gRNA_fate_df[gRNA_fate_df$keep,]

pdf("dist_dist.pdf")
hist(all_dist,breaks = 30)
hist(apply(prop_diff[,c("cDC1_diff","cDC2_diff","pDC_diff")],1,function(x){sqrt(sum(x*x))}),breaks = 30)
dev.off()



library(ggpubr)
cond = c("WT.cDC1", "WT.cDC1","Cas9.cDC1", "Cas9.cDC1","WT.cDC2", "WT.cDC2",  "Cas9.cDC2", "Cas9.cDC2", 
         "WT.pDC", "WT.pDC","Cas9.pDC", "Cas9.pDC")
pdf("sel_genes_fate_vsWT.pdf",width = 10,height = 5)

for (i in gRNA_fate_df$gene_name){
  ix = 1
  pp=list()
  tmp = prop_diff[g_na==i,]
  
  tmp1 = tmp[,c("cDC1_WT",   "cDC2_WT",    "pDC_WT")]
  colnames(tmp1) = c("cDC1","cDC2","pDC")
  tmp1$condition = "WT"
  tmp1$gRNA_id = rownames(tmp1)
  tmp2 = tmp[,c("cDC1_Cas9",   "cDC2_Cas9",    "pDC_Cas9")]
  colnames(tmp2) = c("cDC1","cDC2","pDC")
  tmp2$condition = "Cas9"
  tmp2$gRNA_id = rownames(tmp2)
  
  tmp_df = rbind(tmp1,tmp2)
  row_ix = c()
  for (it in 1:nrow(tmp1)){
    row_ix = c(row_ix,it)
    row_ix = c(row_ix,nrow(tmp1)+it)
  }
  tmp_df = tmp_df[row_ix,]
  rownames(tmp_df) = 1:nrow(tmp_df)
  
  ppp=ggtern(data=tmp_df,aes(cDC1,pDC,cDC2)) + 
    geom_point(aes(col=condition),alpha=0.8,size=2)+
    geom_path(aes(shape=gRNA_id,alpha=0.7),show.legend = F,arrow = arrow(length=unit(0.20,"cm"), type = "closed"))+
    theme_bw()+
    ggtitle(i)
  print(ppp)
  #pp[[ix]] = ggtern(data=tmp,aes(cDC1_Cas9,cDC2_Cas9,pDC_Cas9)) + 
  #  geom_point(alpha=0.8,size=1)+
  #  stat_density_tern(data=prop_WT,aes(cDC1_WT,cDC2_WT,pDC_WT))
  #ix = ix+1
  for (t in 1:nrow(tmp)){
    tmp_df = data.frame(condition=cond,norm_counts=(unlist(tmp[t,1:12])))
    tmp_df$cell_type = as.factor(unlist(lapply(strsplit(as.character(tmp_df$condition),split="[.]"),function(x){x[2]})))
    tmp_df$grp = as.factor(unlist(lapply(strsplit(as.character(tmp_df$condition),split="[.]"),function(x){x[1]})))
    tmp_df$condition = factor(tmp_df$condition, levels = c("WT.cDC1", "Cas9.cDC1", "WT.cDC2","Cas9.cDC2","WT.pDC","Cas9.pDC"))
    pp[[ix]] = ggplot(data=tmp_df,aes(x=cell_type,y=norm_counts,col=grp))+geom_boxplot()+ylim(0,max(tmp_df$norm_counts)+20)+
      stat_compare_means(aes(group = grp), method = "t.test",label = "p.signif")+
      theme_bw()+ggtitle(rownames(tmp[t,]))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ix = ix+1
  }
  print(ggarrange(plotlist=pp,nrow=1,ncol=length(pp),common.legend = TRUE))
}
dev.off()





pdf("sel_genes_fate_vsWT_arrow_top.pdf",width = 10,height = 5)

for (i in sel_genes[sel_genes %in% g_na]){
  ix = 1
  pp=list()
  tmp = prop_diff[g_na==i,]
  
  tmp1 = tmp[,c("cDC1_WT",   "cDC2_WT",    "pDC_WT")]
  colnames(tmp1) = c("cDC1","cDC2","pDC")
  tmp1$condition = "WT"
  tmp1$gRNA_id = rownames(tmp1)
  tmp2 = tmp[,c("cDC1_Cas9",   "cDC2_Cas9",    "pDC_Cas9")]
  colnames(tmp2) = c("cDC1","cDC2","pDC")
  tmp2$condition = "Cas9"
  tmp2$gRNA_id = rownames(tmp2)
  
  tmp_df = rbind(tmp1,tmp2)
  row_ix = c()
  for (it in 1:nrow(tmp1)){
    row_ix = c(row_ix,it)
    row_ix = c(row_ix,nrow(tmp1)+it)
  }
  tmp_df = tmp_df[row_ix,]
  rownames(tmp_df) = 1:nrow(tmp_df)
  
  ppp=ggtern(data=tmp_df,aes(cDC1,pDC,cDC2)) + 
    geom_point(aes(col=condition),alpha=0.8,size=2)+
    geom_path(aes(shape=gRNA_id,alpha=0.7),show.legend = F,arrow = arrow(length=unit(0.20,"cm"), type = "closed"))+
    theme_bw()+
    ggtitle(i)
  print(ppp)
  #pp[[ix]] = ggtern(data=tmp,aes(cDC1_Cas9,cDC2_Cas9,pDC_Cas9)) + 
  #  geom_point(alpha=0.8,size=1)+
  #  stat_density_tern(data=prop_WT,aes(cDC1_WT,cDC2_WT,pDC_WT))
  #ix = ix+1
  for (t in 1:nrow(tmp)){
    tmp_df = data.frame(condition=cond,norm_counts=log2(unlist(tmp[t,1:12])))
    tmp_df$cell_type = as.factor(unlist(lapply(strsplit(as.character(tmp_df$condition),split="[.]"),function(x){x[2]})))
    tmp_df$grp = as.factor(unlist(lapply(strsplit(as.character(tmp_df$condition),split="[.]"),function(x){x[1]})))
    tmp_df$condition = factor(tmp_df$condition, levels = c("WT.cDC1", "Cas9.cDC1", "WT.cDC2","Cas9.cDC2","WT.pDC","Cas9.pDC"))
    pp[[ix]] = ggplot(data=tmp_df,aes(x=cell_type,y=norm_counts,col=grp))+geom_boxplot()+ylim(0,max(tmp_df$norm_counts)+2)+
      stat_compare_means(aes(group = grp), method = "t.test",label = "p.signif")+
      theme_bw()+ggtitle(rownames(tmp[t,]))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ix = ix+1
  }
  print(ggarrange(plotlist=pp,nrow=1,ncol=length(pp)))
}
dev.off()



pdf("sel_genes_fate_vsWT_arrow_second.pdf",width = 10,height = 5)

for (i in sel_snd[sel_snd %in% g_na]){
  ix = 1
  pp=list()
  tmp = prop_diff[g_na==i,]
  
  tmp1 = tmp[,c("cDC1_WT",   "cDC2_WT",    "pDC_WT")]
  colnames(tmp1) = c("cDC1","cDC2","pDC")
  tmp1$condition = "WT"
  tmp1$gRNA_id = rownames(tmp1)
  tmp2 = tmp[,c("cDC1_Cas9",   "cDC2_Cas9",    "pDC_Cas9")]
  colnames(tmp2) = c("cDC1","cDC2","pDC")
  tmp2$condition = "Cas9"
  tmp2$gRNA_id = rownames(tmp2)
  
  tmp_df = rbind(tmp1,tmp2)
  row_ix = c()
  for (it in 1:nrow(tmp1)){
    row_ix = c(row_ix,it)
    row_ix = c(row_ix,nrow(tmp1)+it)
  }
  tmp_df = tmp_df[row_ix,]
  rownames(tmp_df) = 1:nrow(tmp_df)
  
  ppp=ggtern(data=tmp_df,aes(cDC1,pDC,cDC2)) + 
    geom_point(aes(col=condition),alpha=0.8,size=2)+
    geom_path(aes(shape=gRNA_id,alpha=0.7),show.legend = F,arrow = arrow(length=unit(0.20,"cm"), type = "closed"))+
    theme_bw()+
    ggtitle(i)
  print(ppp)
  #pp[[ix]] = ggtern(data=tmp,aes(cDC1_Cas9,cDC2_Cas9,pDC_Cas9)) + 
  #  geom_point(alpha=0.8,size=1)+
  #  stat_density_tern(data=prop_WT,aes(cDC1_WT,cDC2_WT,pDC_WT))
  #ix = ix+1
  for (t in 1:nrow(tmp)){
    tmp_df = data.frame(condition=cond,norm_counts=log2(unlist(tmp[t,1:12])))
    tmp_df$cell_type = as.factor(unlist(lapply(strsplit(as.character(tmp_df$condition),split="[.]"),function(x){x[2]})))
    tmp_df$grp = as.factor(unlist(lapply(strsplit(as.character(tmp_df$condition),split="[.]"),function(x){x[1]})))
    tmp_df$condition = factor(tmp_df$condition, levels = c("WT.cDC1", "Cas9.cDC1", "WT.cDC2","Cas9.cDC2","WT.pDC","Cas9.pDC"))
    pp[[ix]] = ggplot(data=tmp_df,aes(x=cell_type,y=norm_counts,col=grp))+geom_boxplot()+ylim(0,max(tmp_df$norm_counts)+2)+
      stat_compare_means(aes(group = grp), method = "t.test",label = "p.signif")+
      theme_bw()+ggtitle(rownames(tmp[t,]))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ix = ix+1
  }
  print(ggarrange(plotlist=pp,nrow=1,ncol=length(pp)))
}
dev.off()


pdf("test_p_val.pdf")
ggplot(data=tmp_df,aes(x=cell_type,y=norm_counts,col=grp))+geom_boxplot()+ylim(0,max(tmp_df$norm_counts)+2)+
  stat_compare_means(aes(group = grp), method = "t.test",label = "p.signif")+
  theme_bw()+ggtitle(rownames(tmp[t,]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



###### first batch
library(readxl)
X181123_Shalin_Screen_Hits <- read_excel("~/Dropbox/research/sis_seq/validation/CRISPR_screen/181123 Shalin Screen Hits Order - R oligo plate.xlsx")
X181123_Shalin_Screen_Hits = X181123_Shalin_Screen_Hits[-(1:24),]
first_batch_genes = unique(X181123_Shalin_Screen_Hits$Name)


pdf("sel_genes_fate_vsWT_arrow_first_batch.pdf",width = 10,height = 5)

for (i in first_batch_genes[first_batch_genes %in% g_na]){
  ix = 1
  pp=list()
  tmp = prop_diff[g_na==i,]
  
  tmp1 = tmp[,c("cDC1_WT",   "cDC2_WT",    "pDC_WT")]
  colnames(tmp1) = c("cDC1","cDC2","pDC")
  tmp1$condition = "WT"
  tmp1$gRNA_id = rownames(tmp1)
  tmp1$gRNA_txt = ""
  tmp2 = tmp[,c("cDC1_Cas9",   "cDC2_Cas9",    "pDC_Cas9")]
  colnames(tmp2) = c("cDC1","cDC2","pDC")
  tmp2$condition = "Cas9"
  tmp2$gRNA_id = rownames(tmp2)
  tmp2$gRNA_txt = rownames(tmp2)
  
  tmp_df = rbind(tmp1,tmp2)
  row_ix = c()
  for (it in 1:nrow(tmp1)){
    row_ix = c(row_ix,it)
    row_ix = c(row_ix,nrow(tmp1)+it)
  }
  tmp_df = tmp_df[row_ix,]
  rownames(tmp_df) = 1:nrow(tmp_df)
  
  ppp=ggtern(data=tmp_df,aes(cDC1,pDC,cDC2)) + 
    geom_point(aes(col=condition),alpha=0.8,size=2)+
    geom_path(aes(shape=gRNA_id,alpha=0.7),show.legend = F,arrow = arrow(length=unit(0.20,"cm"), type = "closed"))+
    geom_text(aes(label = gRNA_txt), vjust=1,size=3)+
    theme_bw()+
    ggtitle(i)
  print(ppp)
  #pp[[ix]] = ggtern(data=tmp,aes(cDC1_Cas9,cDC2_Cas9,pDC_Cas9)) + 
  #  geom_point(alpha=0.8,size=1)+
  #  stat_density_tern(data=prop_WT,aes(cDC1_WT,cDC2_WT,pDC_WT))
  #ix = ix+1
  for (t in 1:nrow(tmp)){
    tmp_df = data.frame(condition=cond,norm_counts=log2(unlist(tmp[t,1:12])))
    tmp_df$cell_type = as.factor(unlist(lapply(strsplit(as.character(tmp_df$condition),split="[.]"),function(x){x[2]})))
    tmp_df$grp = as.factor(unlist(lapply(strsplit(as.character(tmp_df$condition),split="[.]"),function(x){x[1]})))
    tmp_df$condition = factor(tmp_df$condition, levels = c("WT.cDC1", "Cas9.cDC1", "WT.cDC2","Cas9.cDC2","WT.pDC","Cas9.pDC"))
    pp[[ix]] = ggplot(data=tmp_df,aes(x=cell_type,y=norm_counts,col=grp))+geom_boxplot()+ylim(0,max(tmp_df$norm_counts)+2)+
      stat_compare_means(aes(group = grp), method = "t.test",label = "p.signif")+
      theme_bw()+ggtitle(rownames(tmp[t,]))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ix = ix+1
  }
  print(ggarrange(plotlist=pp,nrow=1,ncol=length(pp)))
}
dev.off()

########

prop_diff_sel = prop_diff[g_na %in% X181123_Shalin_Screen_Hits$Name,]
prop_diff_sel$dist = apply(prop_diff_sel[,c("cDC1_diff", "cDC2_diff", "pDC_diff")],1,function(x){sqrt(sum(x^2))})
prop_diff_sel$max_logFC = apply(log2(prop_diff_sel[,c("cDC1_WT", "cDC2_WT", "pDC_WT")]*prop_diff_sel$WT_sum/(prop_diff_sel[,c("cDC1_Cas9", "cDC2_Cas9", "pDC_Cas9")]*prop_diff_sel$Cas9_sum)),1,function(x){max(abs(x))})


gRNA_lib_sel = gRNA_lib[rownames(prop_diff_sel),]

gRNA_lib_sel$fate_dist = prop_diff_sel$dist
gRNA_lib_sel$max_logFC = prop_diff_sel$max_logFC
gRNA_lib_sel$diff_score = scale(gRNA_lib_sel$fate_dist)+scale(gRNA_lib_sel$max_logFC)

gRNA_lib_sel = gRNA_lib_sel[order(gRNA_lib_sel$Target.Gene.Symbol,gRNA_lib_sel$diff_score,decreasing = TRUE),]
head(gRNA_lib_sel)

priority = 1:4
for (i in 5:nrow(gRNA_lib_sel)){
  if (gRNA_lib_sel$Target.Gene.Symbol[i-1]!=gRNA_lib_sel$Target.Gene.Symbol[i]){
    t = 1
    priority = c(priority,t)
  }else{
    t = t+1
    priority = c(priority,t)
  }
}
gRNA_lib_sel$priority=priority

selected_gRNA_SIS.seq_screen_v1 <- read.csv("~/Dropbox/research/LT03/PromethION/selected_gRNA_SIS-seq_screen_v1.csv", row.names=1, stringsAsFactors=FALSE)

gRNA_lib_sel$gRNA_in_oldbatch="NO"
gRNA_lib_sel$gRNA_in_oldbatch[gRNA_lib_sel$sgRNA.Target.Sequence %in% selected_gRNA_SIS.seq_screen_v1$sgRNA.Target.Sequence]="YES"

write.csv(gRNA_lib_sel,file="selected_gRNA_old_batch_full_table.csv")







sel_second_batch = c("Ctnna3","Drd1","Fcrla","Impdh2","Nrxn1","Siglecg","Tgfb2","Wdr5",
                     "Zfp831","Dip2a","Ifi205","Polm","Tfec","Zbtb14","Zfp455","Zfp652","Zfp831",
                     "Zscan29","Klf12","Ankrd24","Zc3h12a","Naaa","Celsr3","Mbd5","Notch2","Zfp458","Zfy2","Zc3h12a","Mycl")

prop_diff_sel = prop_diff
prop_diff_sel$gene_name = g_na
prop_diff_sel = prop_diff_sel[prop_diff_sel$gene_name %in% c(sel_second_batch,known_marker,sel_mk_genes,"Batf3","Irf8","Id2","Bcor"),]
prop_diff_sel = prop_diff_sel[prop_diff_sel$WT_sum>500 & prop_diff_sel$Cas9_sum>500,]
prop_diff_sel$dist = apply(prop_diff_sel[,c("cDC1_diff", "cDC2_diff", "pDC_diff")],1,function(x){sqrt(sum(x^2))})
prop_diff_sel$max_logFC = apply(log2(prop_diff_sel[,c("cDC1_WT", "cDC2_WT", "pDC_WT")]*prop_diff_sel$WT_sum/(prop_diff_sel[,c("cDC1_Cas9", "cDC2_Cas9", "pDC_Cas9")]*prop_diff_sel$Cas9_sum)),1,function(x){max(abs(x))})

library(dplyr)
tmp = prop_diff_sel[,c("cDC1_diff","cDC2_diff","pDC_diff","gene_name")]

tmp_per_gene = tmp %>% 
  group_by(gene_name) %>% 
  summarise_at(c("cDC1_diff","cDC2_diff","pDC_diff"), function(x){x[which.max( abs(x) )]})
tmp_per_gene = as.data.frame(tmp_per_gene)
rownames(tmp_per_gene) = tmp_per_gene$gene_name
tmp_per_gene = tmp_per_gene[,2:4]
colnames(tmp_per_gene) = c("cDC1","cDC2","pDC")
tmp_per_gene = scale(t(tmp_per_gene))

pheatmap::pheatmap((tmp_per_gene),scale="none",
                   border_color = "grey",
                   color=BlueAndRed(),
                   treeheight_row = 0, treeheight_col = 0,legend=TRUE,fontsize_row=10,
                   fontsize_col=10,
                   filename = "~/Dropbox/research/sis_seq/validation/CRISPR_screen/long_read/heatmap_gene_bias_new.pdf",width = 6,height = 1.5)


prop_diff_sel = prop_diff
prop_diff_sel$gene_name = g_na
prop_diff_sel = prop_diff_sel[prop_diff_sel$gene_name %in% c(sel_second_batch,known_marker,sel_mk_genes,"Batf3","Irf8","Id2","Bcor",unique(c(tmp1_dup, tmp2_dup, tmp3_dup))),]
prop_diff_sel = prop_diff_sel[prop_diff_sel$WT_sum>500 & prop_diff_sel$Cas9_sum>500,]
prop_diff_sel$dist = apply(prop_diff_sel[,c("cDC1_diff", "cDC2_diff", "pDC_diff")],1,function(x){sqrt(sum(x^2))})
prop_diff_sel$max_logFC = apply(log2(prop_diff_sel[,c("cDC1_WT", "cDC2_WT", "pDC_WT")]*prop_diff_sel$WT_sum/(prop_diff_sel[,c("cDC1_Cas9", "cDC2_Cas9", "pDC_Cas9")]*prop_diff_sel$Cas9_sum)),1,function(x){max(abs(x))})

library(dplyr)
tmp = prop_diff_sel[,c("cDC1_diff","cDC2_diff","pDC_diff","gene_name")]

tmp_per_gene = tmp %>% 
  group_by(gene_name) %>% 
  summarise_at(c("cDC1_diff","cDC2_diff","pDC_diff"), function(x){x[which.max( abs(x) )]})
tmp_per_gene = as.data.frame(tmp_per_gene)
rownames(tmp_per_gene) = tmp_per_gene$gene_name
tmp_per_gene = tmp_per_gene[,2:4]
colnames(tmp_per_gene) = c("cDC1","cDC2","pDC")
tmp_per_gene = t(scale(t(tmp_per_gene)))

pheatmap::pheatmap((tmp_per_gene),scale="none",
                   border_color = NA,
                   color=BlueAndRed(),
                   treeheight_row = 0, treeheight_col = 0,legend=TRUE,fontsize_row=4,
                   fontsize_col=10,
                   filename = "~/Dropbox/research/sis_seq/validation/CRISPR_screen/long_read/heatmap_gene_bias_all.pdf",width = 4,height = 9)




gRNA_lib_sel = gRNA_lib[rownames(prop_diff_sel),]

gRNA_lib_sel$fate_dist = prop_diff_sel$dist
gRNA_lib_sel$max_logFC = prop_diff_sel$max_logFC
gRNA_lib_sel$diff_score = scale(gRNA_lib_sel$fate_dist)+scale(gRNA_lib_sel$max_logFC)


gRNA_lib_sel = gRNA_lib_sel[order(gRNA_lib_sel$Target.Gene.Symbol,gRNA_lib_sel$diff_score,decreasing = TRUE),]
head(gRNA_lib_sel)

priority = 1:4
for (i in 5:nrow(gRNA_lib_sel)){
  if (gRNA_lib_sel$Target.Gene.Symbol[i-1]!=gRNA_lib_sel$Target.Gene.Symbol[i]){
    t = 1
    priority = c(priority,t)
  }else{
    t = t+1
    priority = c(priority,t)
  }
}
gRNA_lib_sel$priority=priority

sel_guide = c(1821, 1822, 1823,226, 228,1601, 1603, 1604,933, 
              934, 935,585, 588,2046, 2047,821, 823, 824,1713, 
              1714, 1715, 1716,2321, 2322,1190, 1191, 1192,1869, 1870, 1871, 1872,1049, 1050, 1051, 1052,
              809, 810, 811, 812,913, 914, 915, 916,1843, 1844,2077, 2078, 2080,2321, 2322, 2324,1606, 1607, 1608
              ,429, 431, 432,1337, 1340,1914, 1916,1241, 1243, 1244, 1657, 1658, 1685, 1686,
              569, 570, 571, 572,1989, 1990, 1991, 1992,922,923, 924)

gRNA_lib_sel = gRNA_lib_sel[rownames(gRNA_lib_sel) %in% as.character(sel_guide),]

as.character(sel_guide)[!(as.character(sel_guide) %in% rownames(gRNA_lib_sel))]

write.csv(gRNA_lib_sel,file="selected_gRNA_full_table_no_filter.csv")

write.csv(gRNA_lib_sel[gRNA_lib_sel$priority==1,],file="selected_gRNA_pri_1st.csv")
write.csv(gRNA_lib_sel[gRNA_lib_sel$priority==2,],file="selected_gRNA_pri_2nd.csv")
write.csv(gRNA_lib_sel[gRNA_lib_sel$priority==3,],file="selected_gRNA_pri_3rd.csv")
