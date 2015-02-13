import os, sys, glob

inFilePath = sys.argv[1]
outFilePath = sys.argv[2]

def convertQualityScore(value):
    number = ord(x) - 31
    if number > 73:
        number = 73

    return chr(number)

outFile = open(outFilePath, 'w')
outLines = []

index = 1

for line in file(inFilePath):
    mod = index % 4

    if mod == 0:
        scores = [convertQualityScore(x) for x in line.rstrip()]
        line = "".join(scores) + "\n"

    outLines.append(line)

    if len(outLines) >= 100000:
        outFile.write("".join(outLines))
        outLines = []
        print index

    index += 1

if len(outLines) >= 0:
    outFile.write("".join(outLines))
    print index

outFile.close()
