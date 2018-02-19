library(survival)
setwd('/Users//D/Breast cancer/Okada_20151020')
mutList = read.table('DNB_time5_whole_time5.txt',sep = '\t',header = T,row.names = 1)
mutList = rownames(mutList)

mutList = unlist(lapply(mutList, function(x) unlist(strsplit(x,'[|]'))[1]))
mutList = unique(mutList)

# All patients
AP = read.table('MANIFEST.txt')
AP = AP[,2] # data 2016-01-28 982 samples 
AP = gsub ('-01.maf.txt','',AP)
AP = gsub ('-','.',AP)

mutTable = read.table('MutCount_ BRCA .txt',sep='\t',header = T) # 978 samples

brca.clinical=read.table('gdac.broadinstitute.org_BRCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0/BRCA.clin.merged.picked.txt',header = T,sep = '\t',row.names = 1)
colnames(brca.clinical)=toupper(colnames(brca.clinical))
AP=intersect(colnames(brca.clinical),AP)
brca.clinical=brca.clinical[,AP]


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

data=as.data.frame(cbind(as.numeric(as.character(new.brca.clinical$vital_status)),as.numeric(as.character(new.brca.clinical$survival_days))))
rownames(data)=rownames(new.brca.clinical)
colnames(data)=c('status','days')
mutList = intersect(mutList,unique(mutTable$Gene))
pdf("brca_mutGenes_kmCurve_0,05.pdf")
pval_vec = vector()
for (i in 1:length(mutList)) {
  #i = 3
  gene = mutList[i]
  tmp.table = mutTable[which(mutTable$Gene == gene),]
  tmp.vec = vector()
  MP = apply(tmp.table, 1, function(x) {
    tmp.str = as.character(x[20])
    tmp.str = unlist(strsplit(tmp.str,','))
    #tmp.str = paste(tmp.str,'-01',sep = "")
    tmp.vec = c(tmp.vec,tmp.str)
  })
  MP = unique(gsub(' ','',MP))
  MP = gsub ('-','.',MP)
  NMP = setdiff(AP,MP)
  
  data1 = data[MP,] # mutation sample
  data1$mutation = 1
  data2 = data[NMP,]
  data2$mutation = -1
  data3=rbind(data1,data2)
  
  mfit.byexpr <- survfit(Surv(days, status) ~ mutation, data = data3)
  sdf <- survdiff(Surv(days, status) ~ mutation, data = data3)
  pval <- round(1-pchisq(sdf$chisq, length(sdf$n)-1),3)
  pval_vec = c(pval_vec,pval)
  if (pval <= 0.05) {
    plot(mfit.byexpr,main= paste(gene,'(',length(MP),'Samples)',"in BRCA","pVal=",pval),xlab="Time (Days)", ylab="Overall Survival Proportion",col=c("blue","red"))
    legend("topright",c("Non-mutated","Mutated"),lty = 1:1,col=c("blue","red"))
  }

}
dev.off()

pval.brca.sig = cbind(Gene=mutList,surv.Pval=pval_vec)
write.table(pval.brca.sig,file="brca_mutGene_survival.txt",sep="\t",quote=F,row.names = F)
