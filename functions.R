# parameters----
cms.cols=c("#E79F24","#0071B1","#CA78A6","#009B74")
names(cms.cols) = c("CMS1","CMS2","CMS3","CMS4")
# functions----------------
Combine_Normalize_FF<-function(genes){
  exprs=data.frame(exprs.Marisa[genes,],exprs.TCGA[genes,])
  exprs.nano=data.frame(exprs.nano.AMC90[genes,],exprs.nano.biobank[genes,])
  
  #normalize nanostring-based data together
  library(preprocessCore)
  exprs1=normalize.quantiles(as.matrix(exprs.nano))
  colnames(exprs1)=colnames(exprs.nano)
  rownames(exprs1)=rownames(exprs.nano)
  exprs.nano=exprs1;rm(exprs1)
  
  #normalize other datasets based on nanostring
  for (i in 1:ncol(exprs)) {
    x=normalize.quantiles(as.matrix(data.frame(exprs.nano,exprs[,i])))
    exprs[,i]=x[,ncol(x)]
  }
  
  return(list(exprs=exprs,exprs.nano=exprs.nano))
}
Combine_Normalize_FFPE<-function(genes){
  exprs=data.frame(exprs.Marisa[genes,],exprs.PETACC3[genes,])
  exprs.nano=data.frame(exprs.nano.AMC90[genes,],exprs.nano.biobank[genes,])
  
  #normalize nanostring-based data together
  library(preprocessCore)
  exprs1=normalize.quantiles(as.matrix(exprs.nano))
  colnames(exprs1)=colnames(exprs.nano)
  rownames(exprs1)=rownames(exprs.nano)
  exprs.nano=exprs1;rm(exprs1)
  
  #normalize other datasets based on nanostring
  for (i in 1:ncol(exprs)) {
    x=normalize.quantiles(as.matrix(data.frame(exprs.nano,exprs[,i])))
    exprs[,i]=x[,ncol(x)]
  }
  
  return(list(exprs=exprs,exprs.nano=exprs.nano))
}
bootstraps_function<-function(){
  t(apply(matrix(1:train.repeats, ncol=1),1, 
          function(x) sample(1:nrow(data$col.annot),
                             (nrow(data$col.annot)/n.fold)*(n.fold-1))))  # resampling without replacement
}
alpha_optimizer<-function(alpha=F,w=1){
  #optimize Alpha (the coef for Lasso and Ridge) one time run
  #if optimization needed run: alpha=T
  if(alpha){
    i=1
    exprs.train = cbind(exprs,exprs.nano)
    # select genes by elastic net
    pdf("Alpha.pdf")
    par(mfrow=c(2,2))
    for (j in seq(0,1,0.1)) {
      model=cv.glmnet(t(exprs.train),c(cms,data$col.annot$CMS),
                      family="multinomial", thresh = 1e-07,
                      type.multinomial="grouped", nfolds=10,alpha=j,
                      type.measure = "class",trace.it=0,weights = c(rep(1,ncol(exprs)),rep(w,ncol(exprs.nano))))
      plot(model,main=paste0(j))
      print(j)
    }
    dev.off()
  }
}
design_function<-function(data){
  design = cbind(data$col.annot$CMS,data$col.annot$CMS,data$col.annot$CMS,data$col.annot$CMS); colnames(design) = c("CMS1", "CMS2", "CMS3", "CMS4")
  design[design[,1]!="CMS1",1]=0
  design[design[,2]!="CMS2",2]=0
  design[design[,3]!="CMS3",3]=0
  design[design[,4]!="CMS4",4]=0
  design[design!=0]=1
  design = apply(design,2,as.numeric)
  return(design)
}
sample.results_function<-function(data){
  sample.results =  data.frame(cbind(class=rep(0,nrow(data$col.annot)), "predictions"=0,"errors"=0,"error.rate"=0), 
                               row.names=rownames(data$col.annot))
  sample.results$class = data$col.annot$CMS
  sample.results = cbind(sample.results, design)
  sample.results$CMS1=sample.results$CMS2=sample.results$CMS3=sample.results$CMS4=0
  return(sample.results)
}
train_function<-function(i,b,w){
  return(cv.glmnet(t(exprs.train),c(cms,data$col.annot$CMS[bootstraps[i,]]),
                   family="multinomial", thresh = 1e-07,
                   type.multinomial="grouped", nfolds=10,alpha=alpha,
                   type.measure = "class",#trace.it=1,
                   weights = c(rep(b,ncol(exprs)),rep(w,length(bootstraps[i,])))))
}
model_function<-function(i,b,w){
  tryCatch(train_function(i,b,w) ,
           error=function(e) tryCatch(train_function(i,b,w) ,
                                      error=function(e) tryCatch(train_function(i,b,w) ,
                                                                 error=function(e) tryCatch(train_function(i,b,w) ,
                                                                                            error=function(e) tryCatch(train_function(i,b,w),
                                                                                                                       error=function(e) tryCatch(train_function(i,b,w) ,
                                                                                                                                                  error=function(e) tryCatch(train_function(i,b,w) ,
                                                                                                                                                                             error=function(e) tryCatch(train_function(i,b,w) ,
                                                                                                                                                                                                        error=function(e) model))))))))
}
result_function=function(){
  conf=data.frame(confusion.glmnet(model,newx=t(exprs.test),newy=data$col.annot$CMS[test.samples], family="multinomial"))
  rownames(conf)=paste0(conf[,1],conf[,2])
  correct=conf[c("CMS1CMS1","CMS2CMS2","CMS3CMS3","CMS4CMS4"),"Freq"]
  correct[is.na(correct)]=0
  return(sum(correct)/sum(conf$Freq))
}
update_sample.results<-function(){
  prediction=predict(model, newx=t(exprs.test), interval ="prediction",type="class")
  # update sample.results
  for(j in 1:ncol(exprs.test)){
    sample = colnames(exprs.test)[j]
    if(prediction[j]!=data$col.annot$CMS[test.samples[j]]){ # false prediction
      sample.results[sample,"errors"] = as.numeric(sample.results[sample,"errors"]) + 1
    }
    sample.results[sample,prediction[j]] = sample.results[sample,prediction[j]] + 1
    sample.results[sample,"predictions"] = sample.results[sample,"predictions"] + 1
    sample.results[sample,"error.rate"] = sample.results[sample,"errors"]/sample.results[sample,"predictions"]
  }
  return(sample.results)
}
