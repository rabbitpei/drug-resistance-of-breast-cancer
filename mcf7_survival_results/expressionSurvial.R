#######################################################################
########################### brca survival #############################
#######################################################################

library(survival)
setwd('/Users//jinzeng/Documents/Projects//MCF-7/Results-20161229-Survival/')

target=read.table('T5-12_VS_T1-4_Deseq.txt',sep = '\t',header = T,row.names = 1)
target=read.table('T9-12_VS_T5-8_Deseq.txt',sep = '\t',header = T,row.names = 1)
target=target[which(abs(target$log2FoldChange)>1 & target$padj<0.01),]
target=rownames(target)

brca.rsem=read.table('gdac.broadinstitute.org_BRCA.mRNAseq_Preprocess.Level_3.2016012800.0.0/BRCA.uncv2.mRNAseq_RSEM_all.txt',header = T,sep = '\t',row.names = 1)
tumour=grep("*.01$",colnames(brca.rsem),value=T) 

brca.rsem.tumor=brca.rsem[,colnames(brca.rsem)%in%tumour] # 1093 samples
### remove samples' '.01'
colfile=character()
for (i in 1:ncol(brca.rsem.tumor)){
  tmp=paste(strsplit(colnames(brca.rsem.tumor)[i],split = '[.]')[[1]][1],strsplit(colnames(brca.rsem.tumor)[i],split = '[.]')[[1]][2],strsplit(colnames(brca.rsem.tumor)[i],split = '[.]')[[1]][3],sep = '.')
  colfile=c(colfile,tmp)
}

colnames(brca.rsem.tumor)=colfile

brca.clinical=read.table('gdac.broadinstitute.org_BRCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0/BRCA.clin.merged.picked.txt',header = T,sep = '\t',row.names = 1)
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
for (i in 1:nrow(brca.clinical)){
  if (brca.clinical[i,3]==1) {
    tmp.data=cbind(brca.clinical[i,],brca.clinical[i,4])
    colnames(tmp.data)[ncol(tmp.data)]='survival_days'
    new.brca.clinical=rbind(new.brca.clinical,tmp.data)
  }
  if (brca.clinical[i,3]==0){
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
  data=as.data.frame(cbind(as.numeric(as.character(new.brca.clinical$vital_status)),as.numeric(as.character(new.brca.clinical$survival_days)),brca.rsem.tumor.target[,i]))
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
#write.table(cox.up,file = 'cox_regression.up.txt',quote = F,sep = '\t')

###### plot survival curves
#dir.create("brca_Survival")
#setwd("brca_Survival")
pdf("brca_diffGenesW9-12vsW5-8_kmCurve.pdf")
pval_vec = vector()

#par(mfrow=c(5,3))
for (i in rownames(cox.brca.sig)){
  #i = "A2M"
  data=as.data.frame(cbind(as.numeric(as.character(new.brca.clinical$vital_status)),
                           as.numeric(as.character(new.brca.clinical$survival_days)),
                           brca.rsem.tumor.target[,i]))
  colnames(data)=c('status','days','gene')
  
  data.median=median(data$gene)
  if (data.median != 0) {
    data1=data[which(data$gene>data.median),]
    data1$gene=1
    data2=data[which(data$gene<=data.median),]
    data2$gene=-1
    data3=rbind(data1,data2)
    data3=na.omit(data3)
    mfit.byexpr <- survfit(Surv(days, status) ~ gene, data = data3)
    sdf <- survdiff(Surv(days, status) ~ gene, data = data3)
    pval <- round(1-pchisq(sdf$chisq, length(sdf$n)-1),3)
    
    if (pval <= 0.05) {
      plot(mfit.byexpr,main= paste(i,"in BRCA","pVal=",pval),xlab="Time (Days)", ylab="Overall Survival Proportion",col=c("blue","red"))
      legend("topright",c("low expressed","high expressed"),lty = 1:1,col=c("blue","red"))
    }
  } 
  else {
    pval = 'NaN'
  }
  pval_vec = c(pval_vec,pval)
}
dev.off()

cox.brca.sig = cbind(cox.brca.sig,surv.pval=pval_vec)
write.table(cox.brca.sig,file="brca_diffGenesW9-12vsW5-8_survival.txt",sep="\t",quote=F)
