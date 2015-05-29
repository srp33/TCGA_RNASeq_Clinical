library(pROC)

inFilePath = commandArgs()[7]
actualColumnName = commandArgs()[8]
probabilitiesColumnName = commandArgs()[9]
outFilePath = commandArgs()[10]
main = commandArgs()[11]

data = read.table(inFilePath, sep="\t", stringsAsFactors=F, header=TRUE, row.names=NULL, check.names=F)

actual = as.factor(data[,actualColumnName])
probabilities = as.numeric(data[,probabilitiesColumnName])

pdf(outFilePath)
par(mar=c(4.5, 4.7, 0.0, 0.5), lwd=4)
  
roc_result = roc(actual ~ probabilities, ci=TRUE, plot=TRUE, print.auc=FALSE)
lowerBoundAuc = format(roc_result$ci[1], digits=3)
midAuc = format(roc_result$ci[2], digits=3)
upperBoundAuc = format(roc_result$ci[3], digits=3)
  
ci(roc_result)
sens.ci <- ci.se(roc_result)
plot(sens.ci, type="shape", col="gray95")
plot(sens.ci, type="bars")
plot(roc_result, add=TRUE)
  
text(0.5, 0.00, labels=paste("AUC: ", midAuc, " (", lowerBoundAuc, "-", upperBoundAuc, ")", sep=""))
title(main)
  
par(mar=c(5.1, 4.1, 2.1, 2.1))
graphics.off()
  
print(c(lowerBoundAuc, midAuc, upperBoundAuc))
