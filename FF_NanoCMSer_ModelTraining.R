# load and normalize training data with selected genes-------
source('functions.R')
load("FF.Training.Data.RData")

genes=read.table("FF.Selected.Genes.txt",header = T,stringsAsFactors = F)
genes=genes$Var1
cms=c(cms.Marisa,cms.TCGA)
cms.nano=c(cms.nano.AMC90,cms.nano.biobank)

new_exp=Combine_Normalize_FF(genes)
exprs=new_exp$exprs
exprs.nano=new_exp$exprs.nano

data=list()
data$col.annot=data.frame(row.names = names(cms.nano),CMS=cms.nano)

# train model------------
library(glmnet)
exprs.train = cbind(exprs,exprs.nano)
alpha_optimizer(alpha=T,w=1.5)
alpha=0.7
W=1.5 #Target domain weight

model=cv.glmnet(t(exprs.train),c(cms,data$col.annot$CMS),
                family="multinomial", thresh = 1e-07,
                type.multinomial="grouped", nfolds=10,alpha=alpha,
                type.measure = "class",trace.it=0,
                weights = c(rep(1,ncol(exprs)),rep(W,ncol(exprs.nano))))

save.image("FF_NanoCMSer.RData")


