librarsuppressPackageStartupMessages(library(ROCR))

predictionsFilePath = commandArgs()[7]
outFilePath = commandArgs()[8]
outPdfFilePath = commandArgs()[9]

plotROC = function(actual, probabilities, plotCI=FALSE)
{
  # bottom, left, top, right
  par(mar=c(4.5, 4.7, 0.0, 0.5),lwd=4)

  library(pROC)
  roc_result = roc(actual ~ probabilities, ci=TRUE, plot=TRUE, print.auc=FALSE)
  lowerBoundAuc = format(roc_result$ci[1], digits=3)
  midAuc = format(roc_result$ci[2], digits=3)
  upperBoundAuc = format(roc_result$ci[3], digits=3)

  if (plotCI)
  {
    ci(roc_result)
    sens.ci <- ci.se(roc_result)
    plot(sens.ci, type="shape", col="gray95")
    plot(sens.ci, type="bars")
    plot(roc_result, add=TRUE)
  }

  #text(0.5, 0.00, labels=paste("AUC: ", midAuc, " (", lowerBoundAuc, "-", upperBoundAuc, ")", sep=""))

  par(mar=c(5.1, 4.1, 2.1, 2.1))
}
predictionsFilePath12="~/Desktop/Re__BIOINF-2014-2034_-_Major_Revision/Classification_12_LUAD_LUSC_Predictions.txt"
predictionsFilePath20="~/Desktop/Re__BIOINF-2014-2034_-_Major_Revision/Classification_20_LUAD_LUSC_Predictions.txt"
outPdfFilePath12="~/Desktop/Re__BIOINF-2014-2034_-_Major_Revision/Classification_12_LUAD_LUSC_Predictions.pdf"
outPdfFilePath20="~/Desktop/Re__BIOINF-2014-2034_-_Major_Revision/Classification_20_LUAD_LUSC_Predictions.pdf"
outPngFilePath12="~/Desktop/Re__BIOINF-2014-2034_-_Major_Revision/Classification_12_LUAD_LUSC_Predictions.png"
outPngFilePath20="~/Desktop/Re__BIOINF-2014-2034_-_Major_Revision/Classification_20_LUAD_LUSC_Predictions.png"

outFilePath12="~/Desktop/Re__BIOINF-2014-2034_-_Major_Revision/Classification_12_LUAD_LUSC_Predictions_AUC.txt"
outFilePath20="~/Desktop/Re__BIOINF-2014-2034_-_Major_Revision/Classification_20_LUAD_LUSC_Predictions_AUC.txt"
data12 = read.table(predictionsFilePath12, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)
data20 = read.table(predictionsFilePath20, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1)
actual12 = data12$ActualClass
predictions12 = data12$LUAD_Probability
png(outPngFilePath12)
auc = plotROC(actual12, predictions12, TRUE)
graphics.off()
write(auc, outFilePath12)

actual20 = data20$ActualClass
predictions20 = data20$LUAD_Probability
pdf(outPdfFilePath20)
png(outPngFilePath20)
auc = plotROC(actual20, predictions20, TRUE)
graphics.off()
write(auc, outFilePath20)
