plotHist = function(gene, data, classes, outFilePath)
{
#  xlimMin = min(c(min(data1[gene,]), min(data2[gene,]), min(data3[gene,])))
#  xlimMax = max(c(min(data1[gene,]), max(data2[gene,]), max(data3[gene,])))
#  xlim = c(xlimMin, xlimMax)

  xlim = c(min(data), max(data))
  xlab = paste(gene, " expression levels", sep="")

  pdf(outFilePath)
  par(mfrow=c(3,1),lwd=2)
  hist(as.numeric(data[gene,which(classes=="LUAD")]), main="LUAD", breaks=50, xlab=xlab, xlim=xlim, cex.axis=2,cex.lab=2)
  hist(as.numeric(data[gene,which(classes=="LUSC")]), main="LUSC", breaks=50, xlab=xlab, xlim=xlim,cex.axis=2,cex.lab=2)
  hist(as.numeric(data[gene,which(classes=="Discordant LUSC")]), main="LUSC (potentially discordant)", breaks=12, xlab=xlab, xlim=xlim,cex.axis=2,cex.lab=2)
  graphics.off()
}

library(data.table)
library(heatmap3)
library(RColorBrewer)

data = as.data.frame(fread("RSubread_Discordant_DiffExpressedGenes_Data.txt"))
rownames(data) = data[,1]
data = data.matrix(data[,-1])

classes = read.table("RSubread_Discordant_Classes.txt", sep="\t", stringsAsFactors=F, header=F, row.names=1)

data = data[which(apply(data, 1, var) > 0),]

accent = brewer.pal(8, "Accent")
set3 = brewer.pal(12, "Set3")
cols = c(accent[1], set3[12], accent[7])

ColSideColors = as.character(classes[,1])
ColSideColors[ColSideColors=="LUAD"] = cols[1]
ColSideColors[ColSideColors=="LUSC"] = cols[2]
ColSideColors[ColSideColors=="Discordant LUSC"] = cols[3]
par(lwd=4)
if (nrow(data) <= 10)
  for (gene in rownames(data))
    plotHist(gene, data, classes[,1], paste("RSubread_", gene, "_Histogram.pdf", sep=""))

pdf("RSubread_Discordant_Heatmap.pdf")
colnames(data) = rep("", ncol(data))
if (nrow(data) > 20)
  rownames(data) = rep("", nrow(data))
heatmap3(data, Colv=NA, Rowv=TRUE, showRowDendro=T, showColDendro=F, cexRow=3, margins=c(5, 12), ColSideColors=ColSideColors, ColSideLabs="", cex=1.5)
legend("top", legend=c("LUAD", "LUSC", "Discordant LUSC"), col=cols, cex=1.1, lty=1, lwd=4, inset=-0.07, xpd=TRUE, box.lwd=0, box.lty=0, horiz=F)
graphics.off()
