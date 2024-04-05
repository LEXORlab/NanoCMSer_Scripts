# load and normalize training data-------
source('functions.R')
load("FF.Training.Data.RData")
#remove genes with low correlation
bad.genes=c("SGSM3","TTC31","LY6G6F-LY6G6D","EEF1A1",
            "TRPM4","SPIRE2","BMP1","MAPRE2","RCN3",
            "LMNA","CLCN2","PARP14","NFKB2","GGT5",
            "RMC1","TIMM22","LIMA1","HLA-DRB1","SLITRK6")

genes=setdiff(rownames(exprs.nano.biobank),bad.genes) 
genes=intersect(genes, rownames(exprs.Marisa))
cms=c(cms.Marisa,cms.TCGA)
cms.nano=c(cms.nano.AMC90,cms.nano.biobank)

new_exp=Combine_Normalize_FF(genes)
exprs=new_exp$exprs
exprs.nano=new_exp$exprs.nano

data=list()
data$col.annot=data.frame(row.names = names(cms.nano),CMS=cms.nano)

# ------------------Domain adaptation: optimize weight of target domain--------------
# setting hyperparameters and alpha optimization------------
library(glmnet)
train.repeats=100
genelist.length=200
n.fold=5
W=c(1,0,1,1.5) #W=1.5
bootstraps = bootstraps_function()

alpha_optimizer(alpha=F,w=1)
alpha=0.8 #manually change to the best value by generated figure above

# Optimize Target domain weight-------
start.time = proc.time()[3]

Total=list()
C=T
for (w in W) {
  if(w==1&C){b=0;C=F}else{b=1}
  genelists = matrix(nrow=genelist.length , ncol=train.repeats)
  design = design_function(data)
  sample.results=sample.results_function(data)
  results = vector(length=train.repeats)
  
  for(i in seq(train.repeats)){
    exprs.train = cbind(exprs,exprs.nano[,bootstraps[i,]])
    test.samples = setdiff(1:ncol(exprs.nano), bootstraps[i,])
    exprs.test = exprs.nano[,test.samples]

    model=model_function(i,b,w)
    results[i] =result_function()
    probes=predict(model, type = "coef")[[1]]@i[-1]
    genelists[1:length(probes),i] = probes
    
    sample.results=update_sample.results()
    #if(i%%10==0) print(paste(i, sep=""))
    print(paste("i=",i,"   ","w=",w ,"    ","b=",b,sep=""))
  }
  sample.results = sample.results[order(sample.results$error.rate),]
  Total[[paste0(b,"Weight",w)]]$sample.results=sample.results
  Total[[paste0(b,"Weight",w)]]$results=results
  Total[[paste0(b,"Weight",w)]]$genelists=genelists
}
cat(paste("\n\nrun time training ", (proc.time()[3]-start.time)/(60*60) , " hr\n\n"))


# Result and select optimized weight----------
res=matrix(nrow=train.repeats,ncol=length(W))
Av=rep(0,length(W))
for (i in 1:length(W)) {
  res[,i]=Total[[i]]$results
  Av[i]=1-sum(Total[[i]]$sample.results$errors)/sum(Total[[i]]$sample.results$predictions)
}
names(Av)<-names(Total)

#this value always showed the highest accuracy
ww=which(W%in%1.5)
# ----------------------------------------gene selection----------------------------------
#----------reselect1----------
# Re-combine data after gene selection-----------------
ww=which(W%in%1.5)

gene.freqs = table(Total[[ww]]$genelists[!is.na(Total[[ww]]$genelists)])
gene.freqs = gene.freqs[order(gene.freqs, decreasing=T)]
names(gene.freqs)=rownames(exprs)[as.numeric(names(gene.freqs))] 
genes=names(gene.freqs)[gene.freqs>90]

new_exp=Combine_Normalize_FF(genes)
exprs=new_exp$exprs
exprs.nano=new_exp$exprs.nano

# Reselect the genes--------
alpha_optimizer(alpha=F,w=1.5)
alpha=0.8
train.repeats=1000
genelist.length=200
bootstraps = bootstraps_function()

start.time = proc.time()[3]

genelists = matrix(nrow=genelist.length , ncol=train.repeats)
design = design_function(data)
sample.results=sample.results_function(data)
results = vector(length=train.repeats)

for(i in seq(train.repeats)){
  exprs.train = cbind(exprs,exprs.nano[,bootstraps[i,]])
  test.samples = setdiff(1:ncol(exprs.nano), bootstraps[i,])
  exprs.test = exprs.nano[,test.samples]
  
  model=model_function(i,b,w)
  results[i] =result_function()
  probes=predict(model, type = "coef")[[1]]@i[-1]
  genelists[1:length(probes),i] = probes
  
  sample.results=update_sample.results()
  print(i)
}
sample.results = sample.results[order(sample.results$error.rate),]

cat(paste("\n\nrun time training ", (proc.time()[3]-start.time)/(60*60) , " hr\n\n"))

# accuracy
1-sum(sample.results$errors)/sum(sample.results$predictions)
1-sum(sample.results$error.rate)/sum(sample.results$predictions!=0)

#----------reselect2----------
# Re-combine data after gene selection-----------------
gene.freqs = table(genelists[!is.na(genelists)])
gene.freqs = gene.freqs[order(gene.freqs, decreasing=T)]
names(gene.freqs)=rownames(exprs)[as.numeric(names(gene.freqs))] 

genes=names(gene.freqs)[gene.freqs>900]
new_exp=Combine_Normalize_FF(genes)
exprs=new_exp$exprs
exprs.nano=new_exp$exprs.nano

# Reselect the genes--------
alpha_optimizer(alpha=F,w=1.5)
alpha=0.8
train.repeats=1000
genelist.length=200
bootstraps = bootstraps_function()

start.time = proc.time()[3]

genelists = matrix(nrow=genelist.length , ncol=train.repeats)
design = design_function(data)
sample.results=sample.results_function(data)
results = vector(length=train.repeats)

for(i in seq(train.repeats)){
  exprs.train = cbind(exprs,exprs.nano[,bootstraps[i,]])
  test.samples = setdiff(1:ncol(exprs.nano), bootstraps[i,])
  exprs.test = exprs.nano[,test.samples]
  
  model=model_function(i,b,w)
  results[i] =result_function()
  probes=predict(model, type = "coef")[[1]]@i[-1]
  genelists[1:length(probes),i] = probes
  
  sample.results=update_sample.results()
  print(i)
}
sample.results = sample.results[order(sample.results$error.rate),]

cat(paste("\n\nrun time training ", (proc.time()[3]-start.time)/(60*60) , " hr\n\n"))

# accuracy
1-sum(sample.results$errors)/sum(sample.results$predictions)
1-sum(sample.results$error.rate)/sum(sample.results$predictions!=0)

#----------reselect3----------
# Re-combine data after gene selection-----------------
gene.freqs = table(genelists[!is.na(genelists)])
gene.freqs = gene.freqs[order(gene.freqs, decreasing=T)]
names(gene.freqs)=rownames(exprs)[as.numeric(names(gene.freqs))] 

genes=names(gene.freqs)[gene.freqs>990]
new_exp=Combine_Normalize_FF(genes)
exprs=new_exp$exprs
exprs.nano=new_exp$exprs.nano


# Reselect the genes--------
alpha_optimizer(alpha=F,w=1.5)
alpha=0.8
train.repeats=1000
genelist.length=200
bootstraps = bootstraps_function()

start.time = proc.time()[3]

genelists = matrix(nrow=genelist.length , ncol=train.repeats)
design = design_function(data)
sample.results=sample.results_function(data)
results = vector(length=train.repeats)

for(i in seq(train.repeats)){
  exprs.train = cbind(exprs,exprs.nano[,bootstraps[i,]])
  test.samples = setdiff(1:ncol(exprs.nano), bootstraps[i,])
  exprs.test = exprs.nano[,test.samples]
  
  model=model_function(i,b,w)
  results[i] =result_function()
  probes=predict(model, type = "coef")[[1]]@i[-1]
  genelists[1:length(probes),i] = probes
  
  sample.results=update_sample.results()
  print(i)
}
sample.results = sample.results[order(sample.results$error.rate),]

cat(paste("\n\nrun time training ", (proc.time()[3]-start.time)/(60*60) , " hr\n\n"))

# accuracy
1-sum(sample.results$errors)/sum(sample.results$predictions)
1-sum(sample.results$error.rate)/sum(sample.results$predictions!=0)


#---------------------------------------Save selected genes-----------------------
gene.freqs = table(genelists[!is.na(genelists)])
gene.freqs = gene.freqs[order(gene.freqs, decreasing=T)]
names(gene.freqs)=rownames(exprs)[as.numeric(names(gene.freqs))] 
write.table(gene.freqs, "FF.Selected.Genes.txt",row.names = F)

