# load data-------
source('functions.R')
load("NanoString.Testing.Data.RData")

# Classify RNAseq data using SSP---------
# Run following to install the classifier
# devtools::install_github("Sage-Bionetworks/CMSclassifier")
library(CMSclassifier)
exprs1=log2(exprs.CODE+1)
exprs1=exprs1[!is.na(gene.annot[rownames(exprs1),"entrezId"])&
                !duplicated(gene.annot[rownames(exprs1),"entrezId"]),]
rownames(exprs1)<-gene.annot[rownames(exprs1),"entrezId"]
SSP.CODE <- classifyCMS(exprs1,method="SSP")[[3]]$SSP.predictedCMS

exprs1=log2(exprs.MATCH+1)
exprs1=exprs1[!is.na(gene.annot[rownames(exprs1),"entrezId"])&
                !duplicated(gene.annot[rownames(exprs1),"entrezId"]),]
rownames(exprs1)<-gene.annot[rownames(exprs1),"entrezId"]
SSP.MATCH <- classifyCMS(exprs1,method="SSP")[[3]]$SSP.predictedCMS

# Classify NanoString data using NanoCMSer---------
# Run following to install the classifier
# devtools::install_github("LEXORlab/NanoCMSer")
library("NanoCMSer")
NanoCMSer.CODE <- NanoCMSer(data=exprs.nano.CODE, 
                            sample_type="tumorFFPE", 
                            perform_log2=TRUE, 
                            gene_names='symbol',
                            impute=TRUE)$predictedCMS


NanoCMSer.MATCH <- NanoCMSer(data=exprs.nano.MATCH, 
                             sample_type="tumorFF", 
                             perform_log2=TRUE, 
                             gene_names='symbol',
                             impute=TRUE)$predictedCMS

# consensus matrix-------
mat.c=table(SSP.CODE,NanoCMSer.CODE,useNA = "always")
Conc1=t(apply(mat.c, 1, prop.table))


mat.m=table(SSP.MATCH,NanoCMSer.MATCH,useNA = "always")
Conc2=t(apply(mat.m, 1, prop.table))

pdf("Fig2CD.pdf",3,width = 3.5)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 1), c( "white", "coral3"))
Heatmap(Conc1,cluster_rows = F,cluster_columns = F,
        col=col_fun,rect_gp = gpar(col = "floralwhite", lwd = 2),
        border_gp=gpar(col = "black", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", mat.c[i, j]), x, y)
        })
x=table(NanoCMSer.CODE)
pie(x, labels = paste0(x, "(", round(100 * x/ sum(x)), "%)"),
    col=cms.cols,clockwise=T,main="FFPE-NanoCMSer")

x=table(SSP.CODE)
pie(x, labels = paste0(x, "(", round(100 * x/ sum(x)), "%)"),
    col=cms.cols,clockwise=T,main="SSP")

Heatmap(Conc2,cluster_rows = F,cluster_columns = F,
        col=col_fun,rect_gp = gpar(col = "floralwhite", lwd = 2),
        border_gp=gpar(col = "black", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", mat.m[i, j]), x, y)
        })
x=table(NanoCMSer.MATCH)
pie(x, labels = paste0(x, "(", round(100 * x/ sum(x)), "%)"),
    col=cms.cols,clockwise=T,main="FF-NanoCMSer")

x=table(SSP.MATCH)
pie(x, labels = paste0(x, "(", round(100 * x/ sum(x)), "%)"),
    col=cms.cols,clockwise=T,main="SSP")

dev.off()

