import os, sys, glob
import utilities

inFilePath = sys.argv[1]

inFile = open(inFilePath)
numCols = len(inFile.readline().rstrip().split("\t"))
numRows = 1
for line in inFile:
    numRows += 1
inFile.close()

print "Number Rows: %i" % numRows
print "Number Columns: %i" % numCols
