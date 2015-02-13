import os, sys, glob
from utilities import *

inFilePath = sys.argv[1]
searchPattern = sys.argv[2].decode('string-escape')

print searchPattern in readTextFromFile(inFilePath)
