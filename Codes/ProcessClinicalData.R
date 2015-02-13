
##install packackage 'plyr'
#install.packages("plyr")
library(plyr)

dirname='.'
setwd(dirname)
filenames<-system("ls */nationwidechildrens.org_clinical_patient*", intern=T)

##Identifying only unique column names for all the tumor samples available

for(i in 1:length(filenames)){#####iterating through each of the clinical files to create new matrix files with ALL clinical variables
  f<-read.delim(paste(c(dirname,filenames[i]), collapse=''), header=1) ###reading in the filess one at a time
  rownames(f)<-f$bcr_patient_barcode
  f<-f[3:length(f$bcr_patient_barcode),]
  if(i==1){
    data<-f
  }else{
    data<-list(data,f)
    data<-rbind.fill.matrix(data)
  }
}
rownames(data)<-data[,1]
#Now, converting short TCGA ids reported in clinical data to long TCGA ids reported in RNA-Seq dataset using R codes

pancan20_featureCounts_5<-as.matrix(read.table("Datasets/PANCAN20.IlluminaHiSeq_RNASeqV2.tumor_Rsubread_FeatureCounts.txt", header=1, row.names=1, nrows=5,sep='\t', check.names = F)) 
sample_names<-colnames(pancan20_featureCounts_5)
partial_sample_names<-rownames(data)
counter=0##to check how many replacement has been done
for (j in 1:length(partial_sample_names)){
    if(!is.na(pmatch(partial_sample_names[j],sample_names))){
      partial_sample_names[j]<-sample_names[pmatch(partial_sample_names[j],sample_names, duplicates.ok=F)]  
      counter=counter+1
    }
}
#counter##clinical variables available for 6820 of the 7706 tumor samples
rownames(data)<-partial_sample_names
clinical_data<-matrix(NA, nrow=7706,ncol=420) ###instantiating an NA matrix
rownames(clinical_data)<-colnames(pancan20_featureCounts_5)
colnames(clinical_data)<-colnames(data)
for(i in 1:length(rownames(clinical_data))){
  sample_id<-rownames(clinical_data)[i]
  if(sample_id%in%rownames(data)){
      clinical_data[sample_id,]<-data[sample_id,]      
    }
  }
##Writing the clinical data file

write.table(t(clinical_data),file="TCGA20_clinical_data_ordered_all_clinical_variables_samples_as_columns.txt", sep='\t',col.names=NA, quote=F)
