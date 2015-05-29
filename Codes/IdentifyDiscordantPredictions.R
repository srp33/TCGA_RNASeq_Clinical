inFilePath = commandArgs()[7]
actualColumnName = commandArgs()[8]
predictedColumnName = commandArgs()[9]
potentiallyDiscordantFilePath = commandArgs()[10]

data = read.table(inFilePath, sep="\t", stringsAsFactors=F, header=TRUE, row.names=NULL, check.names=F)

incorrect = data[which(data[,actualColumnName]!=data[,predictedColumnName]),]

potentiallyDiscordantSamples = scan("Potentially_Discordant_LUSC_Samples.txt", what=character(), quiet=TRUE)

print("Samples predicted incorrectly:")
print(nrow(incorrect))

print("Samples predicted incorrectly that were identified previously as potentially discordant:")
print(nrow(incorrect[which(incorrect$row.names %in% potentiallyDiscordantSamples),]))
