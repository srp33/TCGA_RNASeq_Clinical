import os, sys, glob
from utilities import *

inFilePattern = sys.argv[1]
outFilePath = sys.argv[2]

inFilePaths = sorted(glob.glob(inFilePattern))
sampleIDs = [os.path.basename(x) for x in inFilePaths]

features = set()
for inFilePath in inFilePaths:
    print "Identifying features in %s" % inFilePath
    for line in file(inFilePath):
        features.add(line.rstrip().split("\t")[0])
features = sorted(list(features))

outData = [[""] + features]
for inFilePath in inFilePaths:
    print "Parsing and saving values for %s" % inFilePath
    sampleID = os.path.basename(inFilePath)

    valueDict = {}
    for line in file(inFilePath):
        lineItems = line.rstrip().split("\t")
        valueDict[lineItems[0]] = lineItems[1]

    values = [valueDict[feature] for feature in features]
    outData.append([sampleID] + values)

print "Transposing and saving to %s" % outFilePath
writeMatrixToFile(transposeMatrix(outData), outFilePath)
