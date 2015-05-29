library(data.table)
library(stringr)
library(heatmap3)
library(caret)
library(pROC)

readData = function(filePath, logTransform=FALSE)
{
  data = fread(filePath)
  
  data = data.frame(data[-nrow(data),])
  rownames(data) = data[,1]
  data = data[,-1]
  data = data.matrix(data)
  
  if (logTransform)
    data = log2(data + 1)
  
  return(data)
}

mergeData = function(data1, data2)
{
  merged = merge(data1, data2, by=0, sort=FALSE)
  rownames(merged) = merged[,1]
  merged = merged[,-1]
}

crossValidate = function(data, outPrefix)
{
  # Remove any genes with no variance
  data = data[which(apply(data, 1, var) > 0),]
  write.table(dim(data), paste(outPrefix, "_Dimensions.txt", sep=""))

  library(doParallel)
  registerDoParallel(cores=12)

  # From http://stackoverflow.com/questions/13403427/fully-reproducible-parallel-models-using-caret
  # Unfortunately, it doesn't seem to ensure that the results are the same for multiple iterations
  set.seed(0)
  seeds <- vector(mode = "list", length = 11) # length is = (n_repeats*nresampling)+1
  for(i in 1:10) seeds[[i]] <- sample.int(n=1000, 3) #(3 is the number of tuning parameter, mtry for rf, here equal to ncol(iris)-2)
    seeds[[11]]<-sample.int(1000, 1)#for the last model

  model <- train(classes~., data=t(data), method = "rf", tuneLength = 3, trControl = trainControl(method = "cv", savePred=T, classProb=T), seeds=seeds)

  tuneResults = model$results[order(model$results$Accuracy, decreasing=TRUE),]
  bestMtry = tuneResults[1,]$mtry

  # Select predictions that coincide with best mtry parameter
  predictions = model$pred[which(model$pred$mtry == bestMtry),]

  # Sort predictions by the original order
  predictions = predictions[order(predictions$rowIndex),]

  rownames(predictions) = gsub("\\.", "-", colnames(data))

  write.table(predictions, paste(outPrefix, "_Predictions.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=T)

  featureImportance <- varImp(model, scale = TRUE)$importance
  featureImportance <- featureImportance[order(featureImportance$Overall, decreasing=TRUE),,drop=FALSE]
  write.table(featureImportance, paste(outPrefix, "_FeatureImportance.txt", sep=""), quote=FALSE, row.names=T, col.names=NA, sep="\t")
}

identifyDiffExpressedGenes = function(data1, data2, n)
{
  data1Mean = apply(data1, 1, mean)
  data2Mean = apply(data2, 1, mean)
  ratios = (data1Mean + 1) / (data2Mean + 1)
  ratios = sort(ratios, decreasing=TRUE)
  genesToPlot = c(names(head(ratios, n=n)), names(tail(ratios, n=n)))

  return(genesToPlot)
}

tcgaLuad = readData("12_LUAD_t.txt", logTransform=TRUE)
tcgaLusc = readData("12_LUSC_t.txt", logTransform=TRUE)
rsubreadLuad = readData("20_LUAD_t.txt")
rsubreadLusc = readData("20_LUSC_t.txt")

# Extract gene symbols from row names
rownames(tcgaLuad) = sapply(rownames(tcgaLuad), function(x) { str_split(x, "\\|")[[1]][1] })
rownames(tcgaLusc) = sapply(rownames(tcgaLusc), function(x) { str_split(x, "\\|")[[1]][1] })

# Find genes that are common across both data sets
commonTcgaGenes = intersect(rownames(tcgaLuad), rownames(tcgaLusc))
commonRsubreadGenes = intersect(rownames(rsubreadLuad), rownames(rsubreadLusc))
commonGenes = intersect(commonTcgaGenes, commonRsubreadGenes)
nonOverlappingGenes = setdiff(commonRsubreadGenes, commonTcgaGenes)

# Find samples that are common across both data sets
commonLuadSamples = intersect(colnames(tcgaLuad), colnames(rsubreadLuad))
commonLuscSamples = intersect(colnames(tcgaLusc), colnames(rsubreadLusc))

# Select common genes, samples of interest
tcgaLuad = tcgaLuad[commonTcgaGenes,commonLuadSamples]
tcgaLusc = tcgaLusc[commonTcgaGenes,commonLuscSamples]
rsubreadLuad = rsubreadLuad[commonRsubreadGenes,commonLuadSamples]
rsubreadLusc = rsubreadLusc[commonRsubreadGenes,commonLuscSamples]

classesLuad = rep("LUAD", ncol(tcgaLuad))
classesLusc = rep("LUSC", ncol(tcgaLusc))
classes = as.factor(c(classesLuad, classesLusc))

tcga = mergeData(tcgaLuad, tcgaLusc)
rsubread = mergeData(rsubreadLuad, rsubreadLusc)

# Remove any genes with no variance
tcga = tcga[which(apply(tcga, 1, var) > 0),]
rsubread = rsubread[which(apply(rsubread, 1, var) > 0),]

crossValidate(tcga, "TCGA_AllGenes")
crossValidate(rsubread, "RSubread_AllGenes")
crossValidate(tcga[commonGenes,], "TCGA_CommonGenes")
crossValidate(rsubread[commonGenes,], "RSubread_CommonGenes")
crossValidate(rsubread[nonOverlappingGenes,], "RSubread_NonOverlappingGenes")

# Identify top differentially expressed genes
tcgaDiffExpressedGenes = identifyDiffExpressedGenes(tcgaLuad, tcgaLusc, 100)
rsubreadNonOverlappingDiffExpressedGenes = identifyDiffExpressedGenes(rsubreadLuad[nonOverlappingGenes,], rsubreadLusc[nonOverlappingGenes,], 100)

# Get potentially discordant samples
luscDiscordantSamples = scan("Potentially_Discordant_LUSC_Samples.txt", what=character(), quiet=TRUE)
luscDiscordantSamples = str_replace_all(luscDiscordantSamples, "\\-", ".")
luscDiscordantSamples = intersect(luscDiscordantSamples, commonLuscSamples)

tcgaLuscDiscordant = tcgaLusc[,luscDiscordantSamples]
tcgaLusc = tcgaLusc[,setdiff(colnames(tcgaLusc), luscDiscordantSamples)]
tcga = mergeData(tcgaLuad, tcgaLusc)
tcga = mergeData(tcga, tcgaLuscDiscordant)

rsubreadLuscDiscordant = rsubreadLusc[,luscDiscordantSamples]
rsubreadLusc = rsubreadLusc[,setdiff(colnames(rsubreadLusc), luscDiscordantSamples)]
rsubread = mergeData(rsubreadLuad, rsubreadLusc)
rsubread = mergeData(rsubread, rsubreadLuscDiscordant)

#discordantDiffExpressedGenes = identifyDiffExpressedGenes(rsubreadLuad[nonOverlappingGenes,], rsubreadLuscDiscordant[nonOverlappingGenes,], 5)
discordantDiffExpressedGenes = c("MIR320A", "MIR1234", "MIR4461", "MIR186")

colnames(rsubread) = str_replace_all(colnames(rsubread), "\\.", "-")
write.table(rsubread[discordantDiffExpressedGenes,], "RSubread_Discordant_DiffExpressedGenes_Data.txt", sep="\t", quote=F, col.names=NA, row.names=T)

classes = c(classesLuad, rep("LUSC", ncol(rsubreadLusc)), rep("Discordant LUSC", ncol(rsubreadLuscDiscordant)))
classes = cbind(colnames(rsubread), classes)
write.table(classes, "RSubread_Discordant_Classes.txt", sep="\t", quote=F, col.names=F, row.names=F)
