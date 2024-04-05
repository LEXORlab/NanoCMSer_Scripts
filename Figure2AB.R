# load classifier---------
load("All.NanoClassifier.Info.RData")
cms.cols=c("#E79F24","#0071B1","#CA78A6","#009B74"); names(cms.cols) = c("CMS1","CMS2","CMS3","CMS4")
# Fig2A------------
pdf("Fig2A.pdf",height = 3,9)
x=coef(model.FFPE)
x=data.frame(as.matrix(x$CMS1),as.matrix(x$CMS2),as.matrix(x$CMS3),as.matrix(x$CMS4))[-1,]
colnames(x)=paste0("CMS",1:4)
S=apply(x, 1, sd)
x=x[S>0,]
breaks=seq(-0.15,0.15,0.001)
col=colorRampPalette(c("dodgerblue4","white","gold3"))(length(breaks))
library(pheatmap)
pheatmap(t(x),breaks=breaks,col=col,cluster_cols =T,cluster_rows = F,border_color = NA)

x=coef(model.FF)
x=data.frame(as.matrix(x$CMS1),as.matrix(x$CMS2),as.matrix(x$CMS3),as.matrix(x$CMS4))[-1,]
colnames(x)=paste0("CMS",1:4)
S=apply(x, 1, sd)
x=x[S>0,]
breaks=seq(-0.15,0.15,0.001)
col=colorRampPalette(c("dodgerblue4","white","gold3"))(length(breaks))
library(pheatmap)
pheatmap(t(x),breaks=breaks,col=col,cluster_cols =T,cluster_rows = F,border_color = NA,legend = F)
dev.off()


# Fig2B---------------------
pdf("Fig2B.pdf",width = 7,6)
par(mfrow=c(3,4),mar=c(3,2,2,1))
library(scales)
x=coef(model.FFPE)
x=data.frame(as.matrix(x$CMS1),as.matrix(x$CMS2),as.matrix(x$CMS3),as.matrix(x$CMS4))[-1,]
colnames(x)=paste0("CMS",1:4)
y=apply(x,1,function(x){names(x)[x==max(x)]})
genes=c(rownames(which(x==max(x[names(y)[y=="CMS1"],]), arr.ind = TRUE)),
        rownames(which(x==max(x[names(y)[y=="CMS2"],]), arr.ind = TRUE)),
        rownames(which(x==max(x[names(y)[y=="CMS3"],]), arr.ind = TRUE)),
        rownames(which(x==max(x[names(y)[y=="CMS4"],]), arr.ind = TRUE)))

for (j in 1:4) {
  plot(density(exprs.FFPE[genes[j],]),type="n",ylim=c(0,1.5),main=genes[j],ylab="",xlab="",yaxs="i")
  for(i in names(cms.cols)){polygon(density(exprs.FFPE[genes[j],which(cms.FFPE==i)]),col=alpha(cms.cols[i],0.5),border=cms.cols[i])}
  polygon(density(exprs.FFPE[genes[j],which(cms.FFPE==paste0("CMS",j))]),col=alpha(cms.cols[paste0("CMS",j)],0.5),border=cms.cols[paste0("CMS",j)])
  lines(x=c(0,20),y=c(0,0))
}

x=coef(model.FF)
x=data.frame(as.matrix(x$CMS1),as.matrix(x$CMS2),as.matrix(x$CMS3),as.matrix(x$CMS4))[-1,]
colnames(x)=paste0("CMS",1:4)
y=apply(x,1,function(x){names(x)[x==max(x)]})
genes=c(rownames(which(x==max(x[names(y)[y=="CMS1"],]), arr.ind = TRUE)),
        rownames(which(x==max(x[names(y)[y=="CMS2"],]), arr.ind = TRUE)),
        rownames(which(x==max(x[names(y)[y=="CMS3"],]), arr.ind = TRUE)),
        rownames(which(x==max(x[names(y)[y=="CMS4"],]), arr.ind = TRUE)))


for (j in 1:4) {
  plot(density(exprs.FF[genes[j],]),type="n",ylim=c(0,1.5),main=genes[j],ylab="",xlab="",yaxs = "i")
  for(i in names(cms.cols)){polygon(density(exprs.FF[genes[j],which(cms.FF==i)]),col=alpha(cms.cols[i],0.5),border=cms.cols[i])}
  polygon(density(exprs.FF[genes[j],which(cms.FF==paste0("CMS",j))]),col=alpha(cms.cols[paste0("CMS",j)],0.5),border=cms.cols[paste0("CMS",j)])
  lines(x=c(0,20),y=c(0,0))
}

dev.off()
