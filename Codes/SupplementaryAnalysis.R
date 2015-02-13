This is the R code for supplementary analysis

```{r echo=FALSE}
#######################comparing gene counts results in our dataset############
setwd("~/Dropbox/TCGA_PANCAN20_manuscript/Datasets/Analysis datasets/")

TCGA_her2<-read.table("GFP18_HER2_TCGA_normalized_genes.txt", sep='\t', header=1, check.names=F)
nrow(TCGA_her2)
TCGA_her2_filtered<-TCGA_her2[!duplicated(TCGA_her2$Gene),]
rownames(TCGA_her2_filtered)<-TCGA_her2_filtered$Gene
TCGA_her2<-subset(TCGA_her2_filtered,select=-Gene)
tcga_her2_normalized<-data.frame(t(TCGA_her2["ERBB2",]))
colnames(tcga_her2_normalized)<-"TCGA_ERBB2"
sd(tcga_her2_normalized$TCGA_ERBB2)/mean(tcga_her2_normalized$TCGA_ERBB2)

rsub_fpkmlog<-read.table("GFP18_HER2_Rsubread_FPKMlog.txt", sep='\t',header=1, row.names=1, check.names = F)
rsub_fpkmlog_her2<-data.frame(t(rsub_fpkmlog["ERBB2",]))
colnames(rsub_fpkmlog_her2)<-"Rsubread_FPKMlog_ERBB2"
sd(rsub_fpkmlog_her2$Rsubread_FPKMlog_ERBB2)/mean(rsub_fpkmlog_her2$Rsubread_FPKMlog_ERBB2)

rsub_fpkm<-read.table("GFP18_HER2_Rsubread_FPKM.txt", sep='\t',header=1, row.names=1, check.names = F)
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
rsub_tpmlog<-log2(apply(rsub_fpkm,2,fpkmToTpm)+1)
rsub_tpmlog_her2<-data.frame(rsub_tpmlog["ERBB2",])
colnames(rsub_tpmlog_her2)<-"Rsubread_TPMlog_ERBB2"
sd(rsub_tpmlog_her2$Rsubread_TPMlog_ERBB2)/mean(rsub_tpmlog_her2$Rsubread_TPMlog_ERBB2)


library(vioplot)
vioplot(log2(tcga_her2_normalized$TCGA_ERBB2[1:12]+1),log2(tcga_her2_normalized$TCGA_ERBB2[13:17]+1),rsub_fpkmlog_her2$Rsubread_FPKMlog_ERBB2[1:12],rsub_fpkmlog_her2$Rsubread_FPKMlog_ERBB2[13:17],rsub_tpmlog_her2$Rsubread_TPMlog_ERBB2[1:12],rsub_tpmlog_her2$Rsubread_TPMlog_ERBB2[13:17],col='grey')

t.test(log2(tcga_her2_normalized$TCGA_ERBB2[1:12]+1),log2(tcga_her2_normalized$TCGA_ERBB2[13:17]+1))
t.test(rsub_fpkmlog_her2$Rsubread_FPKMlog_ERBB2[1:12],rsub_fpkmlog_her2$Rsubread_FPKMlog_ERBB2[13:17])
t.test(rsub_tpmlog_her2$Rsubread_TPMlog_ERBB2[1:12],rsub_tpmlog_her2$Rsubread_TPMlog_ERBB2[13:17])

sd(log2(tcga_her2_normalized$TCGA_ERBB2))/mean(log2(tcga_her2_normalized$TCGA_ERBB2))
sd(rsub_fpkmlog_her2$Rsubread_FPKMlog_ERBB2)/mean(rsub_fpkmlog_her2$Rsubread_FPKMlog_ERBB2)
sd(rsub_tpmlog_her2$Rsubread_TPMlog_ERBB2)/mean(rsub_tpmlog_her2$Rsubread_TPMlog_ERBB2)
```