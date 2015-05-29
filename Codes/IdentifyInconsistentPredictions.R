inFilePath1 = commandArgs()[7]
inFilePath2 = commandArgs()[8]
actualColumnName = commandArgs()[9]
predictedColumnName = commandArgs()[10]

data1 = read.table(inFilePath1, sep="\t", stringsAsFactors=F, header=TRUE, row.names=NULL, check.names=F)
data2 = read.table(inFilePath2, sep="\t", stringsAsFactors=F, header=TRUE, row.names=NULL, check.names=F)

incorrect1 = data1[which(data1[,actualColumnName]!=data1[,predictedColumnName]),]
incorrect2 = data2[which(data2[,actualColumnName]!=data2[,predictedColumnName]),]

print(nrow(incorrect1))
print(nrow(incorrect2))

diff12 = setdiff(incorrect1$row.names, incorrect2$row.names)
diff21 = setdiff(incorrect2$row.names, incorrect1$row.names)
diffs = c(diff12, diff21)

print("Samples predicted inconsistently between two data sets:")
data = merge(data1, data2, by=1)
print(data[which(data$row.names %in% diffs),])
