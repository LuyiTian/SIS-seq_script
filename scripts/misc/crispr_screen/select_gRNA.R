

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

def_gRNA = broadgpp.brie.library.contents[def %in% broadgpp.brie.library.contents$Target.Gene.Symbol,]

may[!(may %in% broadgpp.brie.library.contents$Target.Gene.Symbol)]

may_gRNA = broadgpp.brie.library.contents[broadgpp.brie.library.contents$Target.Gene.Symbol %in% may,]

may_gRNA = may_gRNA[order(may_gRNA$Rule.Set.2.score),]

