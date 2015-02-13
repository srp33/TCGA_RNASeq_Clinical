#############Reading in predicted HER2 pathway activity################################
setwd(".")
rsub_preds<-read.table("rsubread_10_14.txt", sep='\t', header=1, row.names=1)
head(rsub_preds)
tcga_preds<-read.table("Rsem_10_14.txt", sep='\t', header=1, row.names=1)
all_preds<-merge(rsub_preds,tcga_preds,by=0)
rownames(all_preds)<-all_preds$Row.names
all_preds<-subset(all_preds,select=-Row.names)
clinicals<-t(read.delim('TCGA20_clinical_data_ordered_all_clinical_variables_samples_as_columns.txt',sep='\t',header=1, row.names=1,check.names=F))
brca_clinical<-subset(clinicals,clinicals[,'tumor_tissue_site']=='Breast',select=c("bcr_patient_barcode","her2_status_by_ihc"))
head(brca_clinical)

common_all<-merge(all_preds,brca_clinical,by=0)
#head(common_all)
rownames(common_all)<-common_all$Row.names
common_all<-subset(common_all,select=-Row.names)
head(common_all)

ihc_neg<-subset(common_all,common_all$her2_status_by_ihc=="Negative")#,select=c("Rsubread_fpkm","Rsubread_tpm","TCGA_normalized_gene_counts"))
#length(rownames(ihc_neg))#513
ihc_pos<-subset(common_all,common_all$her2_status_by_ihc=="Positive")#,select=c("Rsubread_fpkm","Rsubread_tpm","TCGA_normalized_gene_counts"))
#length(rownames(ihc_pos))#149

##Checking violin plots for predictions. Violin plots are combinations of boxplot and kernel density of the data, in this case, showing the spread and density of the HER2 pathway activity in TCGA BRCA samples. 

##############violin plot: combination of boxplot and kernel density plot of HER2 pathway activity in TCGA BRCA samples####
library(vioplot)
par(mfrow = c(1, 1))

vioplot(ihc_pos$Rsem_log_q_200_f,ihc_neg$Rsem_log_q_200_f,col='grey', names=c("HER2(+)","HER2(-)"))
title("TCGA RSEM log Pipeline")

vioplot(ihc_pos$FPKM_log_q_200_f,ihc_neg$FPKM_log_q_200_f,col='grey', names=c("HER2(+)","HER2(-)"))
title("Rsubread FPKM log Pipeline")

vioplot(ihc_pos$TPM_log_q_100_f,ihc_neg$TPM_log_q_100_f,col='grey', names=c("HER2(+)","HER2(-)"))
title("Rsubread TPM log Pipeline")


##Now, doing a t-test for comparing the Her2(+) and Her2(-) BRCA samples and measuring coefficient of variation for the predictions

## t-tests comparing HER(+) and HER(-) prediction
t.test(ihc_pos$Rsem_log_q_200_f,ihc_neg$Rsem_log_q_200_f)#p-value = 2.009e-05
t.test(ihc_pos$FPKM_log_q_200_f,ihc_neg$FPKM_log_q_200_f)#p-value = 1.493e-10
t.test(ihc_pos$TPM_log_q_100_f,ihc_neg$TPM_log_q_100_f)#p-value = 3.197e-12

###coefficient of variance analysis##
sd(common_all$Rsem_log_q_200_f)/mean(common_all$Rsem_log_q_200_f)##0.6860363
sd(common_all$FPKM_log_q_200_f)/mean(common_all$FPKM_log_q_200_f)##0.2021242
sd(common_all$TPM_log_q_100_f)/mean(common_all$TPM_log_q_100_f)##0.2898588

