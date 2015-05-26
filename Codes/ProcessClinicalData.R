if (!require("plyr")) {
   install.packages("plyr", dependencies = TRUE)
   library(plyr)
   }

data==identifiers=tmp_data=tmp_identifier=NULL
dirname='.'
setwd(dirname)#Set the directory where the clinical data is located for each cancer in separate folder
filenames<-system("ls */nationwidechildrens.org_clinical_patient*", intern=T)
for(i in 1:length(filenames)){#####iterating through each of the clinical files to create new matrix files with ALL clinical variables
  print(i)
  f<-(read.delim(paste(c(dirname,filenames[i]), collapse=''))) ###reading in the filess one at a time
  tmp_data<-f[3:nrow(f),]
  tmp_identifier<-f[1:2,]
  if(i==1){
    data<-tmp_data
    identifier<-tmp_identifier
}else{
    identifier<-list(identifier,tmp_identifier)
    identifier<-rbind.fill.matrix(identifier)
    for(j in 1:ncol(identifier)){
      if(!is.na(identifier[3,j])){
        identifier[1,j]<-identifier[3,j]
        identifier[2,j]<-identifier[4,j]
      }
    }
    identifier<-identifier[1:2,]
    data<-list(data,tmp_data)
    data<-rbind.fill.matrix(data)
    #data<-merge(data,f)
  }
}
rownames(data)<-data[,2]

#Now, converting short TCGA ids reported in clinical data to long TCGA ids reported in RNA-Seq dataset using R codes

sample_names<-rownames(as.matrix(read.table("PANCAN24_CancerType_Samples.txt", row.names=1, sep='\t', check.names = F))) #getting the long TCGA IDs used in RNA-Seq dataset
partial_sample_names<-rownames(data)
counter=0##to check how many replacement has been done
for (j in 1:length(partial_sample_names)){
  if(!is.na(pmatch(partial_sample_names[j],sample_names))){
    partial_sample_names[j]<-sample_names[pmatch(partial_sample_names[j],sample_names, duplicates.ok=F)]  
    counter=counter+1
  }
}

rownames(data)<-partial_sample_names
clinical_data<-matrix(NA, nrow=9264,ncol=548) ###instantiating an NA matrix
rownames(clinical_data)<-sample_names
colnames(clinical_data)<-colnames(data)
for(i in 1:length(rownames(clinical_data))){
  sample_id<-rownames(clinical_data)[i]
  if(sample_id%in%rownames(data)){
    clinical_data[sample_id,]<-data[sample_id,]      
  }
}
clinical_data_identifier<-cbind(t(identifier),t(clinical_data))
write.table(clinical_data_identifier,file="TCGA_clinical_data_ordered_all_clinical_variables_samples_as_columns.txt", sep='\t',col.names=NA, quote=F)

