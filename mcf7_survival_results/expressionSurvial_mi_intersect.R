#######################################################################
########################### brca survival #############################
#######################################################################

library(survival)
setwd('/Users/liurui/Documents/my documents/MATLAB/breast tumor3/whole dataset/mcf7_survival_results')

DNBs=read.table('DNB_time5_whole_time5',sep='\t',header=F)
target=read.table('T5-12_VS_T1-4_Deseq.txt',sep = '\t',header = T,row.names = 1)
#target=read.table('T9-12_VS_T5-8_Deseq.txt',sep = '\t',header = T,row.names = 1)
target=target[row.names(target)%in%as.character(DNBs[,2]),]   #DNB_target
#target=target[which(abs(target$log2FoldChange)>1 & target$padj<0.01),]
target=rownames(target)

brca.rsem=read.table('BRCA.uncv2.mRNAseq_RSEM_all.txt',header = T,sep = '\t',row.names = 1)
tumour=grep("*.01$",colnames(brca.rsem),value=T) 

brca.rsem.tumor=brca.rsem[,colnames(brca.rsem)%in%tumour] # 1093 samples
### remove samples' '.01'
colfile=character()
for (i in 1:ncol(brca.rsem.tumor)){
  tmp=paste(strsplit(colnames(brca.rsem.tumor)[i],split = '[.]')[[1]][1],strsplit(colnames(brca.rsem.tumor)[i],split = '[.]')[[1]][2],strsplit(colnames(brca.rsem.tumor)[i],split = '[.]')[[1]][3],sep = '.')
  colfile=c(colfile,tmp)
}

colnames(brca.rsem.tumor)=colfile

brca.clinical=read.table('BRCA.clin.merged.picked.txt',header = T,sep = '\t',row.names = 1)
colnames(brca.clinical)=toupper(colnames(brca.clinical))

tmp=intersect(colnames(brca.clinical),colnames(brca.rsem.tumor))
### sort the samples name
brca.rsem.tumor=brca.rsem.tumor[,tmp]
brca.clinical=brca.clinical[,tmp]

##### target genes and data for cox regression
genes=gsub("\\|.+", "", rownames(brca.rsem.tumor))
brca.rsem.tumor$genes=genes
index=duplicated(genes)
brca.rsem.tumor=brca.rsem.tumor[!index,]
rownames(brca.rsem.tumor)=brca.rsem.tumor$genes
brca.rsem.tumor=brca.rsem.tumor[,-ncol(brca.rsem.tumor)]
brca.rsem.tumor.target=brca.rsem.tumor[as.character(target),]

brca.rsem.tumor.target=na.omit(brca.rsem.tumor.target)
brca.rsem.tumor.target=as.data.frame(t(brca.rsem.tumor.target))
brca.clinical=as.data.frame(t(brca.clinical))

new.brca.clinical=data.frame()
for (i in 1:nrow(brca.clinical)) {
  if (brca.clinical[i,3]==1) {
    tmp.data=cbind(brca.clinical[i,],brca.clinical[i,4])
    colnames(tmp.data)[ncol(tmp.data)]='survival_days'
    new.brca.clinical=rbind(new.brca.clinical,tmp.data)
  }
  if (brca.clinical[i,3]==0) {
    tmp.data2=cbind(brca.clinical[i,],brca.clinical[i,5])
    colnames(tmp.data2)[ncol(tmp.data2)]='survival_days'
    new.brca.clinical=rbind(new.brca.clinical,tmp.data2)
  }
}


#########################################################################
############# survival ############################


cox.brca=data.frame()

for (i in 1:ncol(brca.rsem.tumor.target)){
  #i = 1
  data=as.data.frame(cbind(as.numeric(as.character(new.brca.clinical$vital_status)),
                           as.numeric(as.character(new.brca.clinical$survival_days)),
                           brca.rsem.tumor.target[,i]))
  #rownames(data)=tumor$Experiment_Name
  colnames(data)=c('event','days','gene')
  
  a=summary(coxph(formula = Surv(days, event) ~ gene,data))
  cox.new=a$coefficients
  #cox.up.95new=a$conf.int
  rownames(cox.new)=colnames(brca.rsem.tumor.target)[i]
  cox.brca=rbind(cox.brca,cox.new)
  #cox.up.95=rbind(cox.up.95,cox.up.95new)
}

colnames(cox.brca)=c('coef','exp(coef)','se(coef)','z','cox.pval')

#cox.brca.sig=subset(cox.brca,p.value<0.05)
cox.brca.sig = cox.brca
PI=rep(0,nrow(brca.rsem.tumor.target))
for (i in 1:nrow(brca.rsem.tumor.target))   {
  for(j in 1:length(cox.brca.sig$coef))  {
    PI[i]<-PI[i]+brca.rsem.tumor.target[i,j]*(cox.brca.sig$coef[j])
  }
}
target.data<-data
target.data<-cbind(PI,target.data)
medPI<-median(PI)
data1<-target.data[which(target.data$PI>medPI),1:3]
data1$cogene<- 1
data2<-target.data[which(target.data$PI<medPI),1:3]
data2$cogene<- -1
data3=rbind(data1,data2)
#write.table(cox.up,file = 'cox_regression.up.txt',quote = F,sep = '\t')

###### plot survival curves
#dir.create("brca_Survival")
#setwd("brca_Survival")
pdf("brca_diffGenesW9-12vsW5-8_kmCurve_mi_intersect.pdf")
pval_vec = vector()
mfit.byexpr<- survfit(Surv(days,event)~cogene,data=data3)
sdf<-survdiff(Surv(days,event)~cogene,data=data3)
pval<-round(1-pchisq(sdf$chisq,length(sdf$n)-1),3)
plot(mfit.byexpr,main=paste("DNB in Breast cancer", 'pVal=',pval),xlab="Time (Days)",
     ylab="Overrall Survival Proportion", col=c("blue","red"))
dev.off()












#par(mfrow=c(5,3))
#input the gene list***********************input***************************
#select.genes<-read.table('first_survival_SigGenes_JZ.txt')
#select.genes<-as.character(select.genes[,1])
#select.genes<-c('KRT15','HLA-DRB1','AHR','SHB','LARP4B','VHL','ZNF486')
#select.genes<-c('HLA-DRB1','IER3','PSMB8')
select.genes<-c('LARP4B','VHL','ZNF486')
select.genes<-c('VHL')
new.brca.clinical$idx=c(1:dim(new.brca.clinical)[1])
data=as.data.frame(cbind(as.numeric(as.character(new.brca.clinical$idx)),as.numeric(as.character(new.brca.clinical$vital_status)),
                         as.numeric(as.character(new.brca.clinical$survival_days))))
colnames(data)=c('idx','status','days')
samples.med<-array(0,length(select.genes))
l<-0
for (i in select.genes){
  #i = "A2M"
  tit<-
    data$tit<-brca.rsem.tumor.target[,i]
  l<-l+1
  colnames(data)[3+l]<-l
  samples.med[l]<-median(data[,3+l])
}

if (samples.med[1] != 0) {
  positive_set<-as.character(data[which(data[,3+l]>samples.med[l]),]$idx)
  negative_set<-as.character(data[which(data[,3+l]<=samples.med[l]),]$idx)
}

l<-1
for (i in select.genes[-1]) {
  #print(dim(data))
  l<-l+1
  if (samples.med[l] != 0) {
    positive_set<-intersect(positive_set,data[which(data[,3+l]>samples.med[l]),]$idx)
    #data1=data[which(data$gene>tmp$gene[100]),]
    
    negative_set<-intersect(negative_set,data[which(data[,3+l]<=samples.med[l]),]$idx)
    #data2=data[which(data$gene<=tmp$gene[850]),]
  }
}
data1=data[data$idx%in%positive_set==T,1:3]
data1$gene=1
data2=data[data$idx%in%negative_set==T,1:3]
data2$gene=-1

data3=rbind(data1,data2)
data3=na.omit(data3)
mfit.byexpr <- survfit(Surv(days, status) ~ gene, data = data3)
sdf <- survdiff(Surv(days, status) ~ gene, data = data3)
pval <- round(1-pchisq(sdf$chisq, length(sdf$n)-1),3)
pdf("brca_intersect.pdf")
if (pval <= 1) {
  plot(mfit.byexpr,main= paste('INTERSECT',"in BRCA","pVal=",pval),xlab="Time (Days)", ylab="Overall Survival Proportion",col=c("blue","red"))
  legend("topright",c("low expressed","high expressed"),lty = 1:1,col=c("blue","red"))
}
dev.off()

# cox.brca.sig = cbind(cox.brca.sig,surv.pval=pval_vec)
# write.table(cox.brca.sig,file="brca_diffGenesW9-12vsW5-8_survival.txt",sep="\t",quote=F)
