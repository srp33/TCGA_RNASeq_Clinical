find_biological_replicate<-function(matrix){
  s=NULL
  samples=colnames(matrix)
  for(i in 1:ncol(matrix)){
    s[i]=paste(strsplit(samples,'-')[[i]][1:3],sep='',collapse = '-')
  }
  sum(duplicated(s))
  s=s[duplicated(s)]
  dupsamples=NULL
  counter=0
  for(i in 1:length(samples)){ 
    tmp=paste(strsplit(samples[i],'-')[[1]][1:3],sep='',collapse="-")
    print(tmp)
    if(tmp%in%s){
      print("biological replicate found!!")
      print(rownames(samples))[i]
      dupsamples=c(dupsamples,samples[i])
      counter=counter+1
    }
  }

  print(paste(counter,"samples are duplicated for biological replicates"))
  return (matrix[,colnames(matrix)%in%dupsamples])
  
}



library(data.table)

samples<-read.table("~/Downloads/GSE62944_TCGA_20_CancerType_Samples.txt",row.names=1)
# tcga20<-data.frame(fread("~/Desktop/PANCAN20.IlluminaHiSeq_RNASeqV2.tumor_Rsubread_FeatureCounts.txt"),row.names=1,check.names = F)
# dim(tcga20)
# filt_20<-find_biological_replicate(rownames(samples),tcga20)
# dim(filt_20)
# tcga20_zero<-apply(tcga20==0,2,sum)
tcga24_tpm<-data.frame(fread("~/Desktop/PANCAN24/PANCAN24.IlluminaHiSeq_RNASeqV2.tumor_Rsubread_TPM.txt"),row.names=1,check.names = F)
dim(tcga20_tpm)

normal_tpm<-data.frame(fread("~/Desktop/PANCAN24/TCGA24.IlluminaHiSeq_RNASeqV2.normal_Rsubread_TPM.txt"),row.names=1,check.names = F)
dim(normal_tpm)
colnames(normal_tpm)
s=NULL
samples=colnames(normal_tpm)
for(i in 1:ncol(normal_tpm)){
  s[i]=paste(strsplit(samples,'-')[[i]][1:3],sep='',collapse = '-')
}
tumor_s=NULL
samples=colnames(tcga24_tpm)
for(i in 1:ncol(tcga24_tpm)){
  tumor_s[i]=paste(strsplit(samples,'-')[[i]][1:3],sep='',collapse = '-')
}



tcga20_tpm<-data.frame(fread("~/Desktop/PANCAN20.IlluminaHiSeq_RNASeqV2.tumor_Rsubread_TPM_10_9.txt"),row.names=1,check.names = F)
dim(tcga20_tpm)
filt_20_tpm<-find_biological_replicate(tcga20_tpm)
dim(filt_20_tpm)
tcga20_zero<-apply(tcga20_tpm==0,2,sum)


rsem<-data.frame(fread("~/Desktop/PANCAN12.IlluminaHiSeq_RNASeqV2.geneExp.tumor_whitelist"),row.names=1,check.names = F)
filt_12<-find_biological_replicate(rsem)
dim(filt_12)##16 samples with 2 replicates
rsem_zero<-apply(rsem==0,2,sum)
filt_12<-filt_12[30:nrow(filt_12),]
biological_rep_12<-subset(filt_12,select=colnames(filt_12)%in%colnames(filt_20_tpm))
rownames_biological_rep_12<-gsub("[|].*","",rownames(biological_rep_12))

biological_rep_12_o<-biological_rep_12[rownames_biological_rep_12%in%rownames(biological_rep_20),order(colnames(biological_rep_12))]
dim(biological_rep_12_o)

#biological_rep_20<-subset(filt_20,select=colnames(filt_20)%in%colnames(filt_12))
#biological_rep_20_o<-biological_rep_20[,order(colnames(biological_rep_20))]
biological_rep_20<-subset(filt_20_tpm,select=colnames(filt_20_tpm)%in%colnames(filt_12))
biological_rep_20_o<-biological_rep_20[rownames(biological_rep_20)%in%rownames_biological_rep_12,order(colnames(biological_rep_20))]
dim(biological_rep_20_o)


total_20<-log2(apply(biological_rep_20_o,2,sum))
total_12<-log2(apply(biological_rep_12_o,2,sum))


#plot(total_12,ylim=c(22.75,26.5),main="PANCAN12 Level 3",ylab="log2(Gene Counts)")
cor_res_12=cor_res_20=NULL
pdf("~/Dropbox/Bioinformatics submission/Resubmission/scatter_plot.pdf")

cor_res_12=cor_res_20=NULL
for(i in 1:13){
  #points((i*2-1):(i*2),total_12[(i*2-1):(i*2)],col=i,lwd = 4,pch=i)
  plot(log2(biological_rep_12_o[,(i*2-1)]+1),log2(biological_rep_12_o[,(i*2)]+1),xlim=c(0,20),ylim=c(0,20),xlab=paste(colnames(biological_rep_12_o)[(i*2-1)],"log2(Normalized gene counts)",sep='\n'),ylab=paste(colnames(biological_rep_12_o)[(i*2)],"log2(Normalized gene counts)",sep=' '))
  c=cor.test(biological_rep_12_o[,(i*2-1)],biological_rep_12_o[,(i*2)])#,method="spearman")
  cor_res_12=rbind(cor_res_12,c(colnames(biological_rep_12_o)[(i*2-1)],round(total_12[(i*2-1)],digits = 3),colnames(biological_rep_12_o)[(i*2)],round(total_12[(i*2)],digits = 3),round(c$estimate,digits = 3)))
  #title(paste(paste(strsplit(colnames(biological_rep_12_o)[(i*2)],"-")[[1]][1:3],sep="",collapse = "-")," \nrho=",round(c$estimate,digits = 3),sep=""))
  title(paste("TCGA Level 3 \nPearson's correlation=",round(c$estimate,digits = 2),sep=""))
  print(i)
  print(colnames(biological_rep_12_o)[i])
  plot(log2(biological_rep_20_o[,(i*2-1)]+1),log2(biological_rep_20_o[,(i*2)]+1),xlim=c(0,20),ylim=c(0,20),xlab=paste(colnames(biological_rep_20_o)[(i*2-1)],"log2(TPM)",sep='\n'),ylab=paste(colnames(biological_rep_20_o)[(i*2)],"log2(TPM)",sep=''))
  c<-cor.test(biological_rep_20_o[,(i*2-1)],biological_rep_20_o[,(i*2)])#,method="spearman")
  cor_res_20=rbind(cor_res_20,c(colnames(biological_rep_20_o)[(i*2-1)],round(total_20[(i*2-1)],digits = 3),colnames(biological_rep_20_o)[(i*2)],round(total_20[(i*2)],digits=3),round(c$estimate,digits = 3)))
  #title(paste(paste(strsplit(colnames(biological_rep_20_o)[(i*2)],"-")[[1]][1:3],sep="",collapse = "-")," \nrho=",round(c$estimate,digits = 3),sep=""))
  title(paste("Rsubread TPM \nPearson's correlation=",round(c$estimate,digits = 2),sep=""))
  
}
colnames(cor_res_12)=c("Replicate_1","log2 Level 3 gene counts","Replicate_2","log2 Level 3 gene counts","Pearson's correlation between replicates(Level 3)")
colnames(cor_res_20)=c("Replicate_1","log2 Rsubread gene counts","Replicate_2","log2 Rsubread gene counts","Pearson's correlation between replicates(Rsubread)")

par( mfrow = c(2, 1 ) ,lwd=4)
hist(as.numeric(cor_res_12[,5]),main = "TCGA Level 3 Two Replicates\n Each for 13 Samples",xlab = "Pearson's Correlation ", xlim=c(0.88,1),breaks = 5)
abline(v=mean(as.numeric(cor_res_12[,5])),col="red")
abline(v=median(as.numeric(cor_res_12[,5])),col="blue")
hist(as.numeric(cor_res_20[,5]),main = "Rsubread Replicates Two Replicates\n Each for 13 Samples",xlab = "Pearson's Correlation ", xlim=c(0.88,1),breaks = 5)
abline(v=mean(as.numeric(cor_res_20[,5])),col="red")
abline(v=median(as.numeric(cor_res_20[,5])),col="blue")

write.table(cbind(cor_res_12,cor_res_20),"~/Desktop/correlations.txt",sep='\t',col.names=NA,quote=1)
#************************************************************************



####
ecdf_all_ex<-apply(log2(biological_rep_12_o[,c("TCGA-50-5066-01A-01R-1628-07","TCGA-50-5066-02A-11R-2090-07")]+1),2,ecdf)
plot(ecdf_all_ex[[1]],xlab="log2 Level 3 reads", ylab = NA,xlim=c(0,20),col="blue",main="TCGA Level 3",ylim=c(0,1),cex.axis=1.5, cex.lab=1.5)
lines(ecdf_all_ex[[2]],xlab=NA, ylab = NA,col="brown")

###using Rsubread pipeline aligned data
# ecdf_all<-apply(rsub_fpkmlog,2,ecdf)
# plot(ecdf_all[[1]],col="blue",main="Rsubread FPKM",ylim=c(0,1),xlim = c(0,20),cex.axis=1.5, cex.lab=1.5,xlab="log2(normalized expression)",ylab="Cumulative proportion")
# for(i in 2:12){lines(ecdf_all[[i]],xlab=NA,ylab = NA,col="blue")}
# for(i in 13:17){lines(ecdf_all[[i]],xlab=NA,ylab = NA,col="brown")}

ecdf_all_ex<-apply(log2(biological_rep_20_o[,c("TCGA-50-5066-01A-01R-1628-07","TCGA-50-5066-02A-11R-2090-07")]+1),2,ecdf)
plot(ecdf_all_ex[[1]],xlab="log2TPM reads", ylab = NA,xlim=c(0,20),col="blue",main="Rsubread",ylim=c(0,1),cex.axis=1.5, cex.lab=1.5,)
lines(ecdf_all_ex[[2]],xlab=NA, ylab = NA,col="brown")




############zero
setwd("~/Dropbox/TCGA_RNASeq_Clinical/Analysis_datasets/")
rsem_her2_expected_counts<-read.table("GFP18_HER2_TCGA_Pipeline_Expected_Gene_Counts.txt", sep='\t', header=1, row.names=1, check.names=F) 
feature<-read.table("GFP18_HER2_Rsubread_geneCounts.txt", sep='\t',header=1, row.names=1, check.names = F) 
TCGA_her2<-read.table("GFP18_HER2_TCGA_Pipeline_Normalized_Genes_Results.txt", sep='\t', header=1, check.names=F) 
rsub_tpm<-log2(read.table("GFP18_HER2_Rsubread_TPM.txt", sep='\t',header=1, row.names=1, check.names = F)+1) 
TCGA_her2_filtered<-TCGA_her2[!duplicated(TCGA_her2$Gene),]
rownames(TCGA_her2_filtered)<-TCGA_her2_filtered$Gene
TCGA_her2<-subset(TCGA_her2_filtered,select=-Gene)
TCGA_her2_log2<-log2(subset(TCGA_her2_filtered,select=-Gene)+1)

com_genes_TCGA<-TCGA_her2[rownames(TCGA_her2)%in%rownames(rsub_tpm),]
com_genes_TCGA<-com_genes_TCGA[order(rownames(com_genes_TCGA)),]
com_genes_tpm<-rsub_tpm[rownames(rsub_tpm)%in%rownames(com_genes_TCGA),]
com_genes_tpm<-com_genes_tpm[order(rownames(com_genes_tpm)),]
zero_genes_rsem<-com_genes_TCGA[apply(com_genes_TCGA[,1:12]==0,1,mean)!=0,1:12]#atleast one zero in 12 GFP replicates
sum_zero_genes_rsem<-mean(apply(zero_genes_rsem==0,1,sum))##average of how many replicates have same zero expression


nrow(zero_genes_rsem)
nrow(zero_genes_rsem)/nrow(com_genes_TCGA)
mean(apply(zero_genes_rsem,1,mean))#228.859 if TCGA counts are used
zero_genes_f<-com_genes_tpm[apply(com_genes_tpm[,1:12]==0,1,mean)!=0,1:12]##at least one zero in 12 GFP replicates
sum_zero_genes_feature<-mean(apply(zero_genes_f==0,1,sum))##average of how many replicates have same zero expression

nrow(zero_genes_f)
nrow(zero_genes_f)/nrow(com_genes_tpm)
mean(apply(zero_genes_f,1,mean))#0.55 if Rsubread counts are used.
par( mfrow = c(1, 2 ) ,lwd=4)
hist(apply(zero_genes_rsem[,1:12]==0,2,sum),ylim=c(0,10),xlim=c(3900,6500),xlab="Total number of zero expressed \ngene counts per sample",main="TCGA Level 3",breaks=12)
abline(v=median(apply(zero_genes_rsem[,1:12]==0,2,sum)),col="red",lty=2)

hist(apply(zero_genes_f[,1:12]==0,2,sum),ylim=c(0,10),xlim=c(3900,6500),xlab="Total number of zero expressed \ngene counts per sample",main="Rsubread TPM",breaks=2)
abline(v=median(apply(zero_genes_f[,1:12]==0,2,sum)),col="red",lty=2)
pro_t<-zero_sum/nrow(rsem_f)
prop<-cbind(pro_t,pro_r)
colnames(prop)<-c("TCGA Level 3","Rsubread TPM")
write.table(prop,"~/Dropbox/Bioinformatics submission/Resubmission/zero_prop.txt",sep='\t',col.names = NA,quote=F)

print(paste("Total number of nonzero rsubread but zero expressing TCGA genes:",nrow(com_genes_tpm[apply(com_genes_TCGA[,1:12]==0,1,mean)==1&apply(com_genes_tpm[,1:12]>0,1,mean)==1,]),sep=" "))
print(paste("Total number of 1-100 reads in rsubread but zero expressing TCGA genes:",nrow(com_genes_tpm[apply(com_genes_TCGA[,1:12]==0,1,mean)==1&apply(com_genes_tpm[,1:12]>0,1,mean)==1&apply(com_genes_tpm[,1:12]<=100,1,mean)==1,]),sep=" "))
print(paste("Total number of 101-1000 rsubreads but zero expressing TCGA genes:",nrow(com_genes_tpm[apply(com_genes_TCGA[,1:12]==0,1,mean)==1&apply(com_genes_tpm[,1:12]>100,1,mean)==1&apply(com_genes_tpm[,1:12]<=1000,1,mean)==1,]),sep=" "))
print(paste("Total number of 1001-10000 rsubreads but zero expressing TCGA genes:",nrow(com_genes_tpm[apply(com_genes_TCGA[,1:12]==0,1,mean)==1&apply(com_genes_tpm[,1:12]>1000,1,mean)==1&apply(com_genes_tpm[,1:12]<=10000,1,mean)==1,]),sep=" "))
print(paste("Total number of 10000+ rsubreads but zero expressing TCGA genes:",nrow(com_genes_tpm[apply(com_genes_TCGA[,1:12]==0,1,mean)==1&apply(com_genes_tpm[,1:12]>10000,1,mean)==1,]),sep=" "))

print(paste("Total number of nonzero TCGA reads but zero expressing Rsubread genes:",nrow(com_genes_TCGA[apply(com_genes_tpm[,1:12]==0,1,mean)==1&apply(com_genes_TCGA[,1:12]>0,1,mean)==1,]),sep=" "))
print(paste("Total number of 1-100 reads in TCGA but zero expressing Rsubread genes:",nrow(com_genes_TCGA[apply(com_genes_tpm[,1:12]==0,1,mean)==1&apply(com_genes_TCGA[,1:12]>0,1,mean)==1&apply(com_genes_TCGA[,1:12]<=100,1,mean)==1,]),sep=" "))
print(paste("Total number of 101-1000 TCGA but zero expressing Rsubread genes:",nrow(com_genes_TCGA[apply(com_genes_tpm[,1:12]==0,1,mean)==1&apply(com_genes_TCGA[,1:12]>100,1,mean)==1&apply(com_genes_TCGA[,1:12]<=1000,1,mean)==1,]),sep=" "))
print(paste("Total number of 1001-10000 TCGA but zero expressing Rsubread genes:",nrow(com_genes_TCGA[apply(com_genes_tpm[,1:12]==0,1,mean)==1&apply(com_genes_TCGA[,1:12]>1000,1,mean)==1&apply(com_genes_TCGA[,1:12]<=10000,1,mean)==1,]),sep=" "))
print(paste("Total number of 10000+ TCGA but zero expressing Rsubread genes:",nrow(com_genes_TCGA[apply(com_genes_tpm[,1:12]==0,1,mean)==1&apply(com_genes_TCGA[,1:12]>10000,1,mean)==1,]),sep=" "))





# feature_f<-feature[rownames(feature)%in%rownames(rsem_her2_expected_counts),]
# rsem_f<-rsem_her2_expected_counts[rownames(rsem_her2_expected_counts)%in%rownames(feature),]
# zero_sum_r<-apply(feature_f==0,2,sum)
# pro_r<-zero_sum_r/nrow(feature_f)
# pro_r<-pro_r[order(names(pro_r))]
# zero_sum_tcga<-apply(rsem_f==0,2,sum)/nrow(rsem_f)
# print(paste("Total number of zero expressing genes in Rsubread data=",nrow(feature_f[(apply(feature_f[,1:17]==0,1,mean)==1),]),sep=" "))
# print(paste("Total number of zero expressing genes in TCGA data=",nrow(rsem_f[(apply(rsem_f[,1:17]==0,1,mean)==1),]),sep=" "))
# # print(paste("Total number of 1-100 expressing genes in Rsubread data=",nrow(feature_f[(apply(feature_f[,1:17]<100&feature_f[,1:17]<0,1,mean)==1),]),sep=" "))
# # print(paste("Total number of 1-100 expressing genes in TCGA data=",nrow(rsem_f[(apply(rsem_f[,1:17]<100&rsem_f[,1:17]!=0,1,mean)==1),]),sep=" "))
# # print(paste("Total number of 1-100 expressing genes in Rsubread data=",nrow(feature_f[(apply(feature_f[,1:17]<100&feature_f[,1:17]<0,1,mean)==1),]),sep=" "))
# # print(paste("Total number of 1-100 expressing genes in TCGA data=",nrow(rsem_f[(apply(rsem_f[,1:17]<100&rsem_f[,1:17]!=0,1,mean)==1),]),sep=" "))
# 
# dim(rsem_f)
# dim(feature_f)
# rsem_f_o<-rsem_f[order(rownames(rsem_f)),]
# feature_f_o<-feature_f[order(rownames(feature_f)),]
# head(rownames(rsem_f_o))
# head(rownames(feature_f_o))
# zero_genes_rsem<-rsem_f_o[apply(rsem_f_o[,1:17]==0,1,mean)==1&apply(feature_f_o[,1:17]==0,1,mean)==1,]#common zero expressing genes
# nrow(zero_genes_rsem)
# nrow(feature_f_o[apply(feature_f_o[,1:17]==0,1,mean)==1&apply(rsem_f_o[,1:17]==0,1,mean)==1,])
# zero_genes_f<-feature_f_o[apply(rsem_f_o[,1:17]==0,1,mean)!=1&apply(feature_f_o[,1:17]==0,1,mean)==1,]##gene that are zero expressing in feature counts but nonzero in TCGA
# nrow(zero_genes_f)
# zero_genes_r<-rsem_f_o[apply(rsem_f_o[,1:17]==0,1,mean)==1&apply(feature_f_o[,1:17]==0,1,mean)!=1,]##gene that are zero expressing in Level 3 but nonzero in feature
# nrow(zero_genes_r)
# -------
# zero_genes_rsem<-rsem_f_o[apply(rsem_f_o[,1:12]==0,1,mean)!=0,1:12]#atleast one zero in 12 GFP replicates
# nrow(zero_genes_rsem)
# nrow(zero_genes_rsem)/nrow(rsem_f_o)
# mean(apply(zero_genes_rsem,1,mean))
# zero_genes_f<-feature_f_o[apply(feature_f_o[,1:12]==0,1,mean)!=0,1:12]##at least one zero in 12 GFP replicates
# nrow(zero_genes_f)
# nrow(zero_genes_f)/nrow(feature_f_o)
# mean(apply(zero_genes_f,1,mean))

print(paste("Total number of nonzero rsubread but zero expressing TCGA genes:",nrow(feature_f_o[apply(rsem_f_o[,1:12]==0,1,mean)==1&apply(feature_f_o[,1:12]>0,1,mean)==1,]),sep=" "))
print(paste("Total number of 1-100 reads in rsubread but zero expressing TCGA genes:",nrow(feature_f_o[apply(rsem_f_o[,1:12]==0,1,mean)==1&apply(feature_f_o[,1:12]>0,1,mean)==1&apply(feature_f_o[,1:12]<=100,1,mean)==1,]),sep=" "))
print(paste("Total number of 101-1000 rsubreads but zero expressing TCGA genes:",nrow(feature_f_o[apply(rsem_f_o[,1:12]==0,1,mean)==1&apply(feature_f_o[,1:12]>100,1,mean)==1&apply(feature_f_o[,1:12]<=1000,1,mean)==1,]),sep=" "))
print(paste("Total number of 1001-10000 rsubreads but zero expressing TCGA genes:",nrow(feature_f_o[apply(rsem_f_o[,1:12]==0,1,mean)==1&apply(feature_f_o[,1:12]>1000,1,mean)==1&apply(feature_f_o[,1:12]<=10000,1,mean)==1,]),sep=" "))
print(paste("Total number of 10000+ rsubreads but zero expressing TCGA genes:",nrow(feature_f_o[apply(rsem_f_o[,1:12]==0,1,mean)==1&apply(feature_f_o[,1:12]>10000,1,mean)==1,]),sep=" "))

print(paste("Total number of nonzero TCGA reads but zero expressing Rsubread genes:",nrow(rsem_f_o[apply(feature_f_o[,1:12]==0,1,mean)==1&apply(rsem_f_o[,1:12]>0,1,mean)==1,]),sep=" "))
print(paste("Total number of 1-100 reads in TCGA but zero expressing Rsubread genes:",nrow(rsem_f_o[apply(feature_f_o[,1:12]==0,1,mean)==1&apply(rsem_f_o[,1:12]>0,1,mean)==1&apply(rsem_f_o[,1:12]<=100,1,mean)==1,]),sep=" "))
print(paste("Total number of 101-1000 TCGA but zero expressing Rsubread genes:",nrow(rsem_f_o[apply(feature_f_o[,1:12]==0,1,mean)==1&apply(rsem_f_o[,1:12]>100,1,mean)==1&apply(rsem_f_o[,1:12]<=1000,1,mean)==1,]),sep=" "))
print(paste("Total number of 1001-10000 TCGA but zero expressing Rsubread genes:",nrow(rsem_f_o[apply(feature_f_o[,1:12]==0,1,mean)==1&apply(rsem_f_o[,1:12]>1000,1,mean)==1&apply(rsem_f_o[,1:12]<=10000,1,mean)==1,]),sep=" "))
print(paste("Total number of 10000+ TCGA but zero expressing Rsubread genes:",nrow(rsem_f_o[apply(feature_f_o[,1:12]==0,1,mean)==1&apply(rsem_f_o[,1:12]>10000,1,mean)==1,]),sep=" "))





#interesting_genes<-feature_f[rownames(feature_f)%in%rownames(zero_genes),]
#zero_genes_feature<-interesting_genes[apply(interesting_genes[,1:17]==0,1,mean)==1,]#
par( mfrow = c(2, 1 ) ,lwd=4)
hist(apply(zero_genes_rsem[,1:12]==0,2,sum),ylim=c(0,10),xlim=c(3900,6500),xlab="Total number of zero expressed \ngenes counts per sample",main="TCGA",breaks=12)
abline(v=median(apply(zero_genes_rsem[,1:12]==0,2,sum)),col="red",lty=2)

hist(apply(zero_genes_f[,1:12]==0,2,sum),ylim=c(0,10),xlim=c(3900,6500),xlab="Total number of zero expressed \ngenes counts per sample",main="Rsubread",breaks=2)
abline(v=median(apply(zero_genes_f[,1:12]==0,2,sum)),col="red",lty=2)
pro_t<-zero_sum/nrow(rsem_f)
prop<-cbind(pro_t,pro_r)
colnames(prop)<-c("TCGA","Rsubread")
write.table(prop,"~/Dropbox/Bioinformatics submission/Resubmission/zero_prop.txt",sep='\t',col.names = NA,quote=F)
#########LUSC but LUAD-like analysis##
class_12<-read.table("~/Desktop/TCGA_RNASeq_Clinical/Analysis_datasets/Classification_12_LUAD_LUSC_Predictions.txt", header=1,row.names=1)
class_20<-read.table("~/Desktop/TCGA_RNASeq_Clinical/Analysis_datasets/Classification_20_LUAD_LUSC_Predictions.txt", header=1,row.names=1)
mismatch12<-class_12[class_12[,1]!=class_12[,2],]
mismatch20<-class_20[class_20[,1]!=class_20[,2],]
mismatches_all<-merge(class_12[rownames(class_12)%in%rownames(mismatch12)|rownames(class_12)%in%rownames(mismatch20),],class_12[rownames(class_12)%in%rownames(mismatch12)|rownames(class_12)%in%rownames(mismatch20),],by=0)
rownames(mismatches_all)<-gsub("01A-.*-07","01",mismatches_all$Row.names)
mismatches_all<-mismatches_all[,2:ncol(mismatches_all)]
colnames(mismatches_all)<-c("ActualClass.TCGA","PredictedClass.TCGA","LUAD_Probability.TCGA", "LUSC_Probability.TCGA", "ActualClass.Rsubread","PredictedClass.Rsubread",  
                           "LUAD_Probability.Rsubread","LUSC_Probability.Rsubread")
lusc_but_luad<-read.table("~/Dropbox/TCGA_RNASeq_Clinical/Analysis_datasets/LUSC_but_LUAD_like.txt",sep='\t', header=1)
discord<-merge(mismatches_all,lusc_but_luad,by.x=0,by.y=1,all.y=T)
mismatches_all[gsub("01A-.*-07","01",mismatches_all$Row.names)%in%lusc_but_luad$sample,]#identifies the missclassified LUSC, but LUAD-like samples identified by 
lusc_but_luad[!lusc_but_luad$sample%in%gsub("01A-.*-07","01",mismatches_all$Row.names),]
lusc_but_luad[lusc_but_luad$sample%in%gsub("01A-.*-07","01",mismatches_all$Row.names),]
lusc_but_luad[lusc_but_luad$sample%in%gsub("01A-.*-07","01",rownames(mismatch20)),]

