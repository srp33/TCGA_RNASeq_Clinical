import os, sys

inFilePath = sys.argv[1]

file, ext = os.path.splitext(inFilePath)

print ext
