import os, sys, glob
from utilities import *

inFilePattern = sys.argv[1]
outFilePath = sys.argv[2]

outFile = open(outFilePath, 'w')

for inFilePath in glob.glob(inFilePattern):
    outFile.write("%s\t%s\n" % (os.path.basename(inFilePath), readScalarFromFile(inFilePath)))

outFile.close()
