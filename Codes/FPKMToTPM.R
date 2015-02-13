##Reading in FPKM data....
pancan20_fpkm<-read.table("~/Dropbox/TCGA_PANCAN20_manuscript/Datasets/PANCAN20.IlluminaHiSeq_RNASeqV2.tumor_Rsubread_RPKM.txt", sep='\t', header=1, row.names=1,check.names = F)
##Creating a function to convert FPKM values to TPM values

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
##Generating TPM values by calling the function and writing the values to a file
pancan20_TPM<-apply(pancan20_fpkm,2,fpkmToTpm)
write.table(pancan20_TPM,"PANCAN20.IlluminaHiSeq_RNASeqV2.tumor_Rsubread_TPM.txt", quote=F,sep='\t', col.names=NA)
write.table(log2(pancan20_TPM+1),"PANCAN20.IlluminaHiSeq_RNASeqV2.tumor_Rsubread_TPMlog.txt", quote=F,sep='\t', col.names=NA)

