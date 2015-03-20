
##Manually drownload Pancan12 RNA_Seq dataset from https://www.synapse.org/#!Synapse:syn1695324 and Rsubread TPM RNA_Seq data from GEO accession number GSM1536837
pan12<-read.table("PANCAN12.IlluminaHiSeq_RNASeqV2.geneExp.tumor_whitelist", header=1,row.names=1)
pan20<-read.table("GSM1536837_TCGA_20.Illumina.tumor_Rsubread_TPM.txt",header=1,row.names=1)

pan12_f<-pan12[rownames(pan12)%in%rownames(pan20),colnames(pan12)%in%colnames(pan20)]
pan20_f<-pan20[rownames(pan20)%in%rownames(pan12),colnames(pan20)%in%colnames(pan12)]



write.table(apply((pan12_f==0),2,sum),"PANCAN12_19583_by_3380_numZeroes.txt",sep='\t',col.names=F,quote=F)
write.table(apply((pan20_f==0),2,sum),"PANCAN20_19583_by_3380_numZeroes.txt",sep='\t',col.names=F,quote=F)
