import os, sys, glob

inFilePath = sys.argv[1]
outDirPath = sys.argv[2]

inFile = open(inFilePath)
sampleIDs = inFile.readline().rstrip().split("\t")[1:]

lineCount = 0

for line in inFile:
    lineItems = line.rstrip().split("\t")
    gene = lineItems.pop(0)

    for sampleID in sampleIDs:
        outFile = open(outDirPath + "/" + sampleID, 'a')
        outFile.write("%s\t%s\n" % (gene, lineItems.pop(0)))
        outFile.close()

    lineCount += 1
    if lineCount % 1000 == 0:
        print lineCount

inFile.close()
