library(pROC)

inFilePath = commandArgs()[7]
actualColumnName = commandArgs()[8]
predColumnName = commandArgs()[9]

data = read.table(inFilePath, sep="\t", stringsAsFactors=F, header=TRUE, row.names=NULL, check.names=F)

actual = data[,actualColumnName]
pred = data[,predColumnName]

accuracy = sum(actual == pred) / nrow(data)

print(accuracy)
