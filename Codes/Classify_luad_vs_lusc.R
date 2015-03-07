library(caret)

outFilePath12 = "Classification_12_LUAD_LUSC_Predictions.txt"
outFilePath20 = "Classification_20_LUAD_LUSC_Predictions.txt"
# Read data from file
setwd("Analysis_datasets")
luad12=read.table("12_LUAD_t.txt", sep="\t", stringsAsFactors=F, header=TRUE, row.names=1, check.names=F)
lusc12=read.table("12_LUSC_t.txt", sep="\t", stringsAsFactors=F, header=TRUE, row.names=1, check.names=F)
lu12=cbind(luad12,lusc12)
luad20=read.table("20_LUAD_t.txt", sep="\t", stringsAsFactors=F, header=TRUE, row.names=1, check.names=F)
lusc20=read.table("20_LUSC_t.txt", sep="\t", stringsAsFactors=F, header=TRUE, row.names=1, check.names=F)
lu20=cbind(luad20,lusc20)


# only keep the same number of the samples in TCGA processed versus Rsubread processed data
lu20_f<-lu20[,colnames(lu20)%in%colnames(lu12)]

# Remove class values from data frame "LGG"==rownames(data)[9752]
classes12 = as.factor(as.character(lu12[nrow(lu12),]))
data12 = t(data.matrix(lu12[-nrow(lu12),]))
classes20 = as.factor(as.character(lu20_f[nrow(lu20_f),]))
data20 = t(data.matrix(lu20_f[-nrow(lu20_f),]))


# Retain features that do not have zero variance
data12 = data12[,which(apply(data12, 2, var) > 0)]
data20 = data20[,which(apply(data20, 2, var) > 0)]


# Set random seed so results are same each time
set.seed(0)

# Build the classification model
mod12 <- train(classes12~., data=data12, method = "rf", tuneLength = 3, trControl = trainControl(method = "cv", savePred=T, classProb=T))
set.seed(0)
mod20 <- train(classes20~., data=data20, method = "rf", tuneLength = 3, trControl = trainControl(method = "cv", savePred=T, classProb=T))


# Determine which mtry parameter value performed best
tuneResults12 = mod12$results[order(mod12$results$Accuracy, decreasing=TRUE),]
bestMtry12 = tuneResults12[1,]$mtry
tuneResults20 = mod20$results[order(mod20$results$Accuracy, decreasing=TRUE),]
bestMtry20 = tuneResults20[1,]$mtry


# Select predictions that coincide with best mtry parameter
predictions12 = mod12$pred[which(mod12$pred$mtry == bestMtry12),]
predictions20 = mod20$pred[which(mod20$pred$mtry == bestMtry20),]
# Sort predictions by the original order
predictions12 = predictions12[order(predictions12$rowIndex),]
predictions20 = predictions20[order(predictions20$rowIndex),]
# Build output matrix
output12 = cbind(rownames(data12), predictions12[,2], predictions12[,1], predictions12[,3:(ncol(predictions12) - 3)])
colnames(output12) = c("SampleID", "ActualClass", "PredictedClass", paste(colnames(output12)[4:ncol(output12)], "Probability", sep="_"))
output20 = cbind(rownames(data20), predictions20[,2], predictions20[,1], predictions20[,3:(ncol(predictions20) - 3)])
colnames(output20) = c("SampleID", "ActualClass", "PredictedClass", paste(colnames(output20)[4:ncol(output20)], "Probability", sep="_"))


# Save predictions to output file
write.table(output12, outFilePath12, sep="\t", col.names=T, row.names=F, quote=F)
write.table(output20, outFilePath20, sep="\t", col.names=T, row.names=F, quote=F)

