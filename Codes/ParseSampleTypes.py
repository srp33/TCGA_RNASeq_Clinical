import os, sys, glob

inFilePath = sys.argv[1]

sampleTypes = set()
for line in file(inFilePath):
    sampleTypes.add(line.rstrip()[13:15])
print sampleTypes
