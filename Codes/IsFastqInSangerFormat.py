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
        break

    index += 1

print minScore < 40
