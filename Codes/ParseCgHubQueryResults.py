import os, sys, glob

inFilePath = sys.argv[1]
sampleFilePath = sys.argv[2]
outDownloadSamplesDirPath = sys.argv[3]
outCancerTypesDirPath = sys.argv[4]

def parseTagValue(lines, key):
    for line in lines:
        line = line.strip()

        if line.startswith("<%s>" % key):
            return line.replace("/", "").replace("<%s>" % key, "")

    return None

def saveOutput(outLines):
    legacyID = parseTagValue(outLines, "legacy_sample_id")

    if sampleFilePath == "" or legacyID in samplesToKeep:
        analysisID = parseTagValue(outLines, "analysis_id")

        if analysisID != None:
            outFilePath = "%s/%s" % (outDownloadSamplesDirPath, legacyID)
            if os.path.exists(outFilePath):
                print "%s already exists" % outFilePath
            outFile = open(outFilePath, 'w')
            outFile.write("%s\n" % analysisID)
            outFile.close()

            cancerType = parseTagValue(outLines, "disease_abbr")
            if cancerType == None:
                print "Cancer type was not specified for %s." % analysisID
                exit(1)
            outFile = open("%s/%s" % (outCancerTypesDirPath, legacyID), 'w')
            outFile.write("%s\n" % cancerType)
            outFile.close()

            return legacyID

    return None

inFileLines = [line for line in file(inFilePath)]

headerLine1 = inFileLines.pop(0)
headerLine2 = inFileLines.pop(0)

if "Query" in inFileLines[0]:
    inFileLines.pop(0)
    inFileLines.pop(0)

footerLine = inFileLines.pop(len(inFileLines)-1)

if sampleFilePath != "":
    samplesToKeep = set([line.rstrip() for line in file(sampleFilePath)])

samplesSaved = set()

outLines = []

for line in inFileLines:
    if "<Result" in line:
        sampleSaved = saveOutput(outLines)
        if sampleSaved != None:
            samplesSaved.add(sampleSaved)

        if len(samplesSaved) % 100 == 0:
            print "Done processing %i samples" % len(samplesSaved)

        outLines = []

    outLines.append(line)

sampleSaved = saveOutput(outLines)
if sampleSaved != None:
    samplesSaved.add(sampleSaved)
print "Done processing %i samples" % len(samplesSaved)
