# load and normalize training data-------
source('functions.R')
load("All.Training.Data.RData")

# ------------CodSet1-----------------
# combine data-----------------
genes=intersect(intersect(rownames(exprs.nano.AMC90),rownames(exprs.Marisa)),intersect(rownames(exprs.TCGA),rownames(exprs.PETACC3)))
exprs=data.frame(exprs.Marisa[genes,],exprs.TCGA[genes,],exprs.PETACC3[genes,])
exprs.nano=data.frame(exprs.nano.AMC90[genes,],exprs.nano.biobank[genes,])
cms=c(cms.Marisa,cms.TCGA,cms.PETACC3)
cms.nano=c(cms.nano.AMC90,cms.nano.biobank)

library(preprocessCore)
exprs1=normalize.quantiles(as.matrix(exprs.nano))
colnames(exprs1)=colnames(exprs.nano)
rownames(exprs1)=rownames(exprs.nano)
exprs.nano=exprs1;rm(exprs1)


for (i in 1:ncol(exprs)) {
  x=normalize.quantiles(as.matrix(data.frame(exprs.nano,exprs[,i])))
  exprs[,i]=x[,ncol(x)]
}
rm(x)

# heatmaps------------
pdf("Fig1B.pdf",height=4,5)
exprs.z=exprs.TCGA[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
  print(i)
}
breaks=seq(-2,2,.1)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
x=pheatmap::pheatmap(exprs.z[,order(cms.TCGA)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                     annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.TCGA),
                     cluster_cols = F,show_colnames = F,show_rownames = F,main = paste0("TCGA, n = ",ncol(exprs.z)))

pheatmap::pheatmap(exprs.z[x$tree_row$order,order(cms.TCGA)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                   annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.TCGA),
                   cluster_cols = F,show_colnames = F,show_rownames = F,cluster_rows = F,main = paste0("TCGA, n = ",ncol(exprs.z)),legend=F,annotation_legend = F)

exprs.z=exprs.Marisa[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
  print(i)
}
breaks=seq(-2,2,.1)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
pheatmap::pheatmap(exprs.z[x$tree_row$order,order(cms.Marisa)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                   annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.Marisa),
                   cluster_cols = F,show_colnames = F,show_rownames = F,cluster_rows = F,main =paste0("GSE39582, n = ",ncol(exprs.z)),legend=F,annotation_legend = F)

exprs.z=exprs.PETACC3[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
  print(i)
}
breaks=seq(-2,2,.1)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
pheatmap::pheatmap(exprs.z[x$tree_row$order,order(cms.PETACC3)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                   annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.PETACC3),
                   cluster_cols = F,show_colnames = F,show_rownames = F,cluster_rows = F,main = paste0("PETACC3, n = ",ncol(exprs.z)),legend=F,annotation_legend = F)


exprs.z=exprs.nano[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
  print(i)
}
breaks=seq(-2,2,.1)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
pheatmap::pheatmap(exprs.z[x$tree_row$order,order(cms.nano)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                   annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.nano#,Tissue=gsub(".*\\.","",colnames(exprs.z))
                   ),
                   cluster_cols = F,show_colnames = F,show_rownames = F,cluster_rows = F,border_color = NA,
                   main = paste0("NanoString, n = ",ncol(exprs.z)),legend=F,annotation_legend = F)
dev.off()

# ------------CodSet2-----------------
# combine data-----------------
load("All.NanoClassifier.Info.RData")
genes=unique(c(rownames(exprs.FF),rownames(exprs.FFPE)))
genes=intersect(intersect(genes,rownames(exprs.Marisa)),intersect(rownames(exprs.TCGA),rownames(exprs.PETACC3)))
exprs=data.frame(exprs.Marisa[genes,],exprs.TCGA[genes,],exprs.PETACC3[genes,])
exprs.nano=data.frame(exprs.nano.AMC90[genes,],exprs.nano.biobank[genes,])
cms=c(cms.Marisa,cms.TCGA,cms.PETACC3)
cms.nano=c(cms.nano.AMC90,cms.nano.biobank)

library(preprocessCore)
exprs1=normalize.quantiles(as.matrix(exprs.nano))
colnames(exprs1)=colnames(exprs.nano)
rownames(exprs1)=rownames(exprs.nano)
exprs.nano=exprs1;rm(exprs1)


for (i in 1:ncol(exprs)) {
  x=normalize.quantiles(as.matrix(data.frame(exprs.nano,exprs[,i])))
  exprs[,i]=x[,ncol(x)]
}
rm(x)

#heatmaps------------
pdf("Fig1C.pdf",height=2,5)
exprs.z=exprs.TCGA[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
  print(i)
}
breaks=seq(-2,2,.1)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
x=pheatmap::pheatmap(exprs.z[,order(cms.TCGA)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                     annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.TCGA),
                     cluster_cols = F,show_colnames = F,show_rownames = F,main = paste0("TCGA, n = ",ncol(exprs.z)))

pheatmap::pheatmap(exprs.z[x$tree_row$order,order(cms.TCGA)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                   annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.TCGA),
                   cluster_cols = F,show_colnames = F,show_rownames = F,cluster_rows = F,main = paste0("TCGA, n = ",ncol(exprs.z)),legend=F,annotation_legend = F)

exprs.z=exprs.Marisa[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
  print(i)
}
breaks=seq(-2,2,.1)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
pheatmap::pheatmap(exprs.z[x$tree_row$order,order(cms.Marisa)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                   annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.Marisa),
                   cluster_cols = F,show_colnames = F,show_rownames = F,cluster_rows = F,main =paste0("GSE39582, n = ",ncol(exprs.z)),legend=F,annotation_legend = F)

exprs.z=exprs.PETACC3[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
  print(i)
}
breaks=seq(-2,2,.1)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
pheatmap::pheatmap(exprs.z[x$tree_row$order,order(cms.PETACC3)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                   annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.PETACC3),
                   cluster_cols = F,show_colnames = F,show_rownames = F,cluster_rows = F,main = paste0("PETACC3, n = ",ncol(exprs.z)),legend=F,annotation_legend = F)


exprs.z=exprs.nano[genes,]
for (i in 1:nrow(exprs.z)) {
  exprs.z[i,]=(exprs.z[i,]-mean(as.numeric(exprs.z[i,])))/sd(exprs.z[i,])
  print(i)
}
breaks=seq(-2,2,.1)
col=colorRampPalette(c("#1483E4","white","#FF8C00"))(length(breaks))
pheatmap::pheatmap(exprs.z[x$tree_row$order,order(cms.nano)],annotation_colors = list(CMS=cms.cols),breaks = breaks,color = col,
                   annotation_col = data.frame(row.names=colnames(exprs.z),CMS=cms.nano#,Tissue=gsub(".*\\.","",colnames(exprs.z))
                   ),
                   cluster_cols = F,show_colnames = F,show_rownames = F,cluster_rows = F,border_color = NA,
                   main = paste0("NanoString, n = ",ncol(exprs.z)),legend=F,annotation_legend = F)
dev.off()
