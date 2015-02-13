import os, sys
import utilities

inFilePath = sys.argv[1]
outFilePath = sys.argv[2]

data = utilities.readMatrixFromFile(inFilePath)

if len(data) > 1 and len(data[0]) == len(data[1]) - 1:
    data[0].insert(0, " ")

utilities.writeMatrixToFile(utilities.transposeMatrix(data), outFilePath)
