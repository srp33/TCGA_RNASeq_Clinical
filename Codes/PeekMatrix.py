import os, sys, glob

inFilePath = sys.argv[1]
numRows = int(sys.argv[2])
numCols = int(sys.argv[3])

lineCount = 0
for line in file(inFilePath):
    lineItems = line.rstrip().split("\t")

    if lineCount <  numRows:
        print lineItems[:numCols]
        lineCount += 1
    else:
        break
