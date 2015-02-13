import os, sys, glob

inFilePath = sys.argv[1]

index = 1

maxScore = -1
minScore = 1000

for line in file(inFilePath):
    if index % 4 == 0:
        scores = [ord(x) for x in line.rstrip()]
        if min(scores) < minScore:
            minScore = min(scores)
        if max(scores) > maxScore:
            maxScore = max(scores)

    if index % 100000 == 0:
        print index

    index += 1

print "Min score: %i" % minScore
print "Max score: %i" % maxScore

# See http://en.wikipedia.org/wiki/FASTQ_format
print "The score range for Sanger format is 33 to 126"
print "The score range for Illumina 1.0 format is 59 to 126"
print "The score range for Illumina 1.3 format is 64 to 126"
