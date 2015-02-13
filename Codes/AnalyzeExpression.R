#######################comparing gene counts results in our dataset############
setwd("Datasets")

rsem_her2_expected_counts<-read.table("HMEC_TCGA_pipeline_GFP_HER2_final.txt", sep='\t', header=1, row.names=1, check.names=F)
rsem_her2<-data.frame(t(rsem_her2_expected_counts["ERBB2",]))
colnames(rsem_her2)<-"TCGA_ERBB2"

feature<-read.table("GFP18_HER2_featureCounts.txt", sep='\t',header=1, row.names=1, check.names = F)
rsub_her2<-data.frame(t(feature["ERBB2",]))
colnames(rsub_her2)<-"Rsubread_ERBB2"

par(mfrow = c(1, 1))
boxplot(log2(rsem_her2$TCGA_ERBB2[1:12]),log2(rsem_her2$TCGA_ERBB2[13:17]),log2(rsub_her2$Rsubread_ERBB2[1:12]),log2(rsub_her2$Rsubread_ERBB2[13:17]),col='grey',xlab='',ylab="log2(gene counts)",main=paste('Comparing TGCA and Rsubread Pipelines','\n', 'in Differentiating HER2 Overexpression from Controls',sep=''))
names=c('TCGA GFP','TCGA HER2','Rsubread GFP', 'Rsubread HER2')
text(seq(1,4,by=1),par("usr")[3] - 1, labels = names, srt = 45, pos = 1, xpd = TRUE)

#Now, let's checks t-test and coefficient of variation for the gene counts.
##Using not normalized data processed by RSEM detected difference in her2 gene count in HER2 overexpressed versus GFP overexpressed samples
##t = -12.1833, df = 4.157, p-value = 0.0002081 but was worse than Rsubread

t.test(log2(rsem_her2$TCGA_ERBB2[1:12]),log2(rsem_her2$TCGA_ERBB2[13:17]))

##Using normalized data processed by Rsubread was good at detecting differences in her2 gene count in HER2 overexpressed versus GFP overexpressed samples
##t = -46.6747, df = 8.35, p-value = 2.152e-11
t.test(log2(rsub_her2$Rsubread_ERBB2[1:12]),log2(rsub_her2$Rsubread_ERBB2[13:17]))

##coefficient of variation is lower for Rsubread processed data than the TCGA MapSplice+RSEM pipeline
mtext("TCGA p-value=0.0002081, Rsubread p-value=2.152e-11",side=3)
sd(rsub_her2$Rsubread_ERBB2)/mean(rsub_her2$Rsubread_ERBB2)  ##1.529863 for our pipeline
sd(rsem_her2$TCGA_ERBB2)/mean(rsem_her2$TCGA_ERBB2) ###1.70127 for TCGA pipeline

