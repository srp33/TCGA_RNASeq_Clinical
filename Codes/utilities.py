import glob, os, posix, sys, math, collections, json, difflib
#import scipy
#from scipy.stats import *
from operator import itemgetter, attrgetter
import itertools
from random import uniform, sample
#import numpy
from collections import defaultdict
#from fisher import *
#from transcendental import stdtr

def printFlush(text, outFilePath=None):
    print text
    sys.stdout.flush()

    if outFilePath != None:
        outFile = open(outFilePath, 'a')
        outFile.write(text + "\n")
        outFile.close()

def printMatrix(data):
    for x in data:
        print x
    print ""

def smartDivide(numerator, denominator):
    if float(denominator) == 0.0:
        return float('nan')

    return float(numerator) / float(denominator)

def getProbes(probeTabFilePath):
    probes = []

    probeTabFile = open(probeTabFilePath)
    headerItems = [x.lower() for x in probeTabFile.readline().rstrip().split("\t")]

    for line in probeTabFile:
        lineItems = line.rstrip().split("\t")
        if headerItems.count("probe set name") > 0:
            probeset = lineItems[headerItems.index("probe set name")]
        else:
            if headerItems.count("probe set id") > 0:
                probeset = lineItems[headerItems.index("probe set id")]
            else:
                print "No probe set name or probe set id column in %s" % probeTabFilePath

        probeX = lineItems[headerItems.index("probe x")]
        probeY = lineItems[headerItems.index("probe y")]
        probe = probeset + "#" + probeX + "_" + probeY
        probes.append(probe)

    return probes

def getProbesetProbesDict(probes):
    probesetProbesDict = {}

    for probe in probes:
        probeset = probe[:probe.find("#")]
        probesetProbesDict[probeset] = probesetProbesDict.setdefault(probeset, []) + [probe]

    return probesetProbesDict

def getPatientIDs(normDirPath, normFileSuffix):
    ids = []

    #print normDirPath + "*" + normFileSuffix
    #sys.exit(0)
    for filePath in glob.glob(normDirPath + "*" + normFileSuffix):
        ids.append(filePath.replace(normDirPath, "").replace(normFileSuffix, ""))

    ids.sort()
    return ids

def readScalarFromFile(filePath):
    return readMatrixFromFile(filePath)[0][0]

def writeScalarToFile(x, filePath):
    outFile = open(filePath, 'w')
    outFile.write(x)
    outFile.close()

def readVectorFromFile(filePath):
    return [line.rstrip() for line in file(filePath)]

def writeVectorToFile(data, filePath):
    outFile = open(filePath, 'w')
    for x in data:
        outFile.write(str(x) + "\n")
    outFile.close()

def readMatrixFromFile(filePath, numLines=None):
    matrix = []
    for line in file(filePath):
        if numLines != None and len(matrix) >= numLines:
            break

        matrix.append(line.rstrip().split("\t"))

        if len(matrix) % 100000 == 0:
            print len(matrix)

    return matrix

def writeMatrixToFile(x, filePath, writeMode='w'):
    outFile = open(filePath, writeMode)
    writeMatrixToOpenFile(x, outFile)
    outFile.close()

def writeMatrixToOpenFile(x, outFile):
    for y in x:
        outFile.write("\t".join([str(z) for z in y]) + "\n")

def appendMatrixToFile(x, filePath):
    writeMatrixToFile(x, filePath, writeMode='a')

def readTextFromFile(filePath):
    text = ""

    for line in file(filePath):
        text += line

    return text

def writeDictToFile(dictionary, filePath):
    writeScalarToFile(json.dumps(dictionary), filePath)

def readDictFromFile(filePath):
    txt = readTextFromFile(filePath)
    dictionary = json.loads(txt)

    dictionary2 = {}

    for key in dictionary:
        value = dictionary[key]

        if isNumeric(key):
            key = int(key)

        dictionary2[key] = value

    return dictionary2

def calculateMean(values):
    if len(values) == 0:
        return float('nan')

    return sum(values) / len(values)

def calculateVarianceMean(values):
    mu = calculateMean(values)
    diffValues = [(x - mu)**2 for x in values]
    return calculateMean(diffValues) / (len(diffValues) - 1)

def calculateWeightedMean(values, weights):
    if len(values) != len(weights):
        print "When calculating a weighted mean, the values must be the same length as the weights."
        raise

def calculateStandardDeviation(values):
    xbar = calculateMean(values)
    residuals = [x - xbar for x in values]
    residualsSquared = [x**2 for x in residuals]
    return math.sqrt(sum(residualsSquared) / (len(values) - 1))

def calculateZscore(x):
    mean = calculateMean(x)
    standardDeviation = calculateStandardDeviation(x)
    return [(y - mean) / standardDeviation for y in x]

def calculateTrimmedMean(values, trimProportion=0.10):
    if values == None or len(values) == 0:
        return None

    values = sorted([float(x) for x in values])

    if len(values) < 3:
        return calculateMean(values)
    elif len(values) == 3:
        return values[1]
    elif len(values) == 4:
        return calculateMean(values[1:3])
    elif len(values) == 5:
        return calculateMean(values[1:4])

    values = scipy.stats.trimboth(values, trimProportion)

    return float(calculateMean(values))

def calculateEuclideanDistance(xList, yList):
    zSum = 0.0

    for i in range(len(xList)):
        x = xList[i]
        y = yList[i]
        z = math.pow(x - y, 2)
        zSum += z

    return math.sqrt(zSum)

def calculateCorrelationCoefficient(xList, yList):
    return numpy.corrcoef(xList, yList)[0,1]

def calculatePearsonCoefficient(xList, yList):
    return stats.pearsonr(xList, yList)[0]

def calculateSpearmanCoefficient(xList, yList):
    return stats.spearmanr(xList, yList)[0]

def calculateTTest(xList, yList):
    xList = numpy.array([x for x in xList if not math.isnan(x)])
    yList = numpy.array([y for y in yList if not math.isnan(y)])

    if len(xList) == 1 and len(yList) > 1:
        return calculateOneSampleTTest(xList[0], yList)
    if len(xList) > 1 and len(yList) == 1:
        return calculateOneSampleTTest(yList[0], xList)

    return ttest_ind(xList, yList, 0)[1]

# From http://stackoverflow.com/questions/10038543/tracking-down-the-assumptions-made-by-scipys-ttest-ind-function
def calculateWelchTTest(pop1, pop2):
    num1 = numpy.array(pop1).shape[0]
    num2 = numpy.array(pop2).shape[0]

    t_stat = (numpy.mean(pop1) - numpy.mean(pop2))/numpy.sqrt( numpy.var(pop1)/num1 + numpy.var(pop2)/num2)
    df = ((numpy.var(pop1)/num1 + numpy.var(pop2)/num2)**(2.0)) / ((numpy.var(pop1)/num1)**(2.0)/(num1-1) + (numpy.var(pop2)/num2) ** (2.0) / (num2-1))

    #one_tailed_p_value = 1.0 - scipy.stats.t.cdf(t_stat,df)
    two_tailed_p_value = 1.0 - (scipy.stats.t.cdf(numpy.abs(t_stat),df) - scipy.stats.t.cdf(-numpy.abs(t_stat), df))

    return two_tailed_p_value

def calculateOneSampleTTest(x, yList):
    return stats.ttest_1samp(yList, x)[1]

def isValueAberrant(x, yList, numStandardDeviations):
    std = calculateStandardDeviation(yList)
    lowerLimit = calculateMean(yList) - float(numStandardDeviations) * std
    upperLimit = calculateMean(yList) + float(numStandardDeviations) * std

    return x < lowerLimit or x > upperLimit

def calculateMedian(values):
  sortedValues = sorted(values)

  if len(sortedValues) % 2 == 1:
      return sortedValues[(len(sortedValues)+1)/2-1]
  else:
      lower = sortedValues[len(sortedValues)/2-1]
      upper = sortedValues[len(sortedValues)/2]
      return (float(lower + upper)) / 2

def calculateFoldChange(values1, values2):
    overallMin = min(min(values1), min(values2))

    values1 = [x - overallMin + 1 for x in values1]
    values2 = [x - overallMin + 1 for x in values2]

    mean1 = calculateMean(values1)
    mean2 = calculateMean(values2)

    return mean1 / mean2

def calculateAbsoluteFoldChange(values1, values2):
    overallMin = min(min(values1), min(values2))

    values1 = [x - overallMin + 1 for x in values1]
    values2 = [x - overallMin + 1 for x in values2]

    mean1 = calculateMean(values1)
    mean2 = calculateMean(values2)

    ratioA = mean1 / mean2
    ratioB = mean2 / mean1

    return min(ratioA, ratioB)

def getNormalizedProbes(normFilePath):
    print "Getting normalized probes"
    return [line.split(" ")[0] for line in file(normFilePath)]

def getKeyProbeDict(filePath, probesToKeep=None, minProbesPerKey=1):
    probesToKeepSet = set(probesToKeep)
    keyProbeDict = {}

    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        key = lineItems[0]

        if len(lineItems) > 1:
            fileProbes = [x for x in lineItems[1].split(",") if x != ""]

            if len(fileProbes) >= 0:
                keyProbeDict[key] = keyProbeDict.setdefault(key, []) + fileProbes

    return keyProbeDict

def getTranscriptProbeDict(filePath, normFilePath):
    normalizedProbes = set(getNormalizedProbes(normFilePath))

    print "Getting transcript-probe dictionary"
    transcriptProbeDict = {}
    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        transcript = lineItems[0]
        probes = lineItems[1].split(",")
        probes = list(set(probes) & normalizedProbes)

        transcriptProbeDict[transcript] = probes

    return transcriptProbeDict

def getPatientsKeyValuesDict(sourceDir, patientIDs, fileSuffix, dataValueIndex, keys=None):
    patientsKeyValuesDict = collections.defaultdict(dict)

    if len(patientIDs) == 0:
        return patientsKeyValuesDict

    keyLineIndicesDict = {}
    lineCount = 0

    for line in file(sourceDir + patientIDs[0] + fileSuffix):
        key = line.rstrip().split("\t")[0]
        keyLineIndicesDict[key] = lineCount

        lineCount += 1
        #if lineCount % 100000 == 0:
        #    print "Parsing file line indices: %i" % lineCount

    #print "Creating key line indices list from dict"
    if keys == None:
        keyLineIndices = [(key, keyLineIndicesDict[key]) for key in keyLineIndicesDict.keys()]
    else:
        keyLineIndices = [(key, keyLineIndicesDict[key]) for key in keys if key in keyLineIndicesDict.keys()]

    #print "Sorting key line indices"
    keyLineIndices.sort(key=itemgetter(1))

    patientFileHandles = {}
    for patientID in patientIDs:
        patientFileHandles[patientID] = open(checkDirPath(sourceDir) + patientID + fileSuffix)

    for patientID in patientIDs:
        #print patientID
        patientFile = open(checkDirPath(sourceDir) + patientID + fileSuffix)

        previousLineIndex = 0
        for keyLineIndex in keyLineIndices:
            for i in range(previousLineIndex, keyLineIndex[1]):
                patientFile.readline()
            previousLineIndex = keyLineIndex[1] + 1

            lineItems = patientFile.readline().rstrip().split("\t")
            patientsKeyValuesDict[patientID][lineItems[0]] = lineItems[dataValueIndex]

        patientFile.close()

    return patientsKeyValuesDict

def getPatientKeyValuesDict(filePath, dataColumnIndex, probes=None):
    probeValues = {}

    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        probe = lineItems[0]
        value = lineItems[dataColumnIndex]

        probeValues[probe] = value

    if not probes:
        return probeValues
    else:
        modProbeValues = {}
        for probe in probes:
            modProbeValues[probe] = probeValues[probe]
        return modProbeValues

def savePatientKeyValuesDict(patientDict, outFilePath):
    outFile = open(outFilePath, 'w')

    for key in sorted(patientDict.keys()):
        outFile.write("%s\t%s\n" % (key, patientDict[key]))

    outFile.close()

def checkDirPath(dirPath):
    if not os.path.exists(dirPath):
        posix.mkdir(dirPath)

    if not dirPath.endswith("/"):
        dirPath = dirPath + "/"

    return dirPath

def lastIndexOf(theList, value):
    return len(theList) - 1 - theList[::-1].index(value)

def getTranscriptGeneDict(filePath):
    transcriptGeneDict = {}

    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        transcript = lineItems[0]

        gene = lineItems[1]
        if len(lineItems) == 3:
            gene = lineItems[2]

        transcriptGeneDict[transcript] = gene

    return transcriptGeneDict

def getGeneTranscriptDict(filePath):
    geneTranscriptDict = {}

    for line in file(filePath):
        lineItems = line.rstrip().split("\t")
        transcript = lineItems[0]

        gene = lineItems[1]
        if len(lineItems) == 3:
            gene = lineItems[2]

        geneTranscriptDict[gene] = geneTranscriptDict.setdefault(gene, []) + [transcript]

    return geneTranscriptDict

def transposeMatrix(x):
    transposed = zip(*x)

    for i in range(len(transposed)):
        transposed[i] = list(transposed[i])

    return transposed

# Copied from: http://code.activestate.com/recipes/491268-ordering-and-ranking-for-lists/
def order(x, NoneIsLast = True, decreasing = False):
    """
    Returns the ordering of the elements of x. The list
    [ x[j] for j in order(x) ] is a sorted version of x.

    Missing values in x are indicated by None. If NoneIsLast is true,
    then missing values are ordered to be at the end.
    Otherwise, they are ordered at the beginning.
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True

    n  = len(x)
    ix = range(n)
    if None not in x:
        ix.sort(reverse = decreasing, key = lambda j : x[j])
    else:
        # Handle None values properly.
        def key(i, x = x):
            elem = x[i]
            # Valid values are True or False only.
            if decreasing == NoneIsLast:
                return not(elem is None), elem
            else:
                return elem is None, elem
        ix = range(n)
        ix.sort(key=key, reverse=decreasing)

    if omitNone:
        n = len(x)
        for i in range(n-1, -1, -1):
            if x[ix[i]] == None:
                n -= 1
        return ix[:n]
    return ix

# Copied from: http://code.activestate.com/recipes/491268-ordering-and-ranking-for-lists/
def rankSmart(x, NoneIsLast=True, decreasing = False, ties = "first"):
    """
    Returns the ranking of the elements of x. The position of the first
    element in the original vector is rank[0] in the sorted vector.

    Missing values are indicated by None.  Calls the order() function.
    Ties are NOT averaged by default. Choices are:
                 "first" "average" "min" "max" "random" "average"
    """
    omitNone = False
    if NoneIsLast == None:
        NoneIsLast = True
        omitNone = True
    O = order(x, NoneIsLast = NoneIsLast, decreasing = decreasing)
    R = O[:]
    n = len(O)
    for i in range(n):
        R[O[i]] = i
    if ties == "first" or ties not in ["first", "average", "min", "max", "random"]:
        return R

    blocks     = []
    isnewblock = True
    newblock   = []
    for i in range(1,n) :
        if x[O[i]] == x[O[i-1]]:
            if i-1 not in newblock:
                newblock.append(i-1)
            newblock.append(i)
        else:
            if len(newblock) > 0:
                blocks.append(newblock)
                newblock = []
    if len(newblock) > 0:
        blocks.append(newblock)

    for i, block  in enumerate(blocks):
        # Don't process blocks of None values.
        if x[O[block[0]]] == None:
            continue
        if ties == "average":
            s = 0.0
            for j in block:
                s += j
            s /= float(len(block))
            for j in block:
                R[O[j]] = s
        elif ties == "min":
            s = min(block)
            for j in block:
                R[O[j]] = s
        elif ties == "max":
            s =max(block)
            for j in block:
                R[O[j]] = s
        elif ties == "random":
            s = sample([O[i] for i in block], len(block))
            for i,j in enumerate(block):
                R[O[j]] = s[i]
        else:
            for i,j in enumerate(block):
                R[O[j]] = j
    if omitNone:
        R = [ R[j] for j in range(n) if x[j] != None]
    return R

# The following function came from http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
def rank2(a):
    n = len(a)
    ivec=rank_simple(a)
    svec=[a[rank] for rank in ivec]
    sumranks = 0
    dupcount = 0
    newarray = [0]*n
    for i in xrange(n):
        sumranks += i
        dupcount += 1
        if i==n-1 or svec[i] != svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in xrange(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0

    return newarray

def globFilesSortedByModTime(pattern):
    def getModifiedTime(filename):
        return os.stat(filename).st_mtime

    return sorted(glob.glob(pattern), key=getModifiedTime)

## From http://stackoverflow.com/questions/34518/natural-sorting-algorithm
def naturalSort(x, reverse=False):
    def natural_key(s):
        return tuple(
            int(''.join(chars)) if isdigit else ''.join(chars)
            for isdigit, chars in itertools.groupby(s, str.isdigit)
        )

    return sorted(x, key=natural_key, reverse=reverse)

def getItemFrequencyMap(x):
    d = defaultdict(int)
    for item in x:
        d[item] += 1

    return d

from math import modf, floor

def quantile(x, q,  qtype = 7, issorted = False):
    """
    Args:
       x - input data
       q - quantile
       qtype - algorithm
       issorted- True if x already sorted.

    Compute quantiles from input array x given q.For median,
    specify q=0.5.

    References:
       http://reference.wolfram.com/mathematica/ref/Quantile.html
       http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile

    Author:
    Ernesto P.Adorio Ph.D.
    UP Extension Program in Pampanga, Clark Field.
    """
    if not issorted:
        y = sorted(x)
    else:
        y = x
    if not (1 <= qtype <= 9):
       return None  # error!

    # Parameters for the Hyndman and Fan algorithm
    abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3

            (0,   0, 0, 1), # California linear interpolation, R type 4
            (0.5, 0, 0, 1), # hydrologists method, R type 5
            (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
            (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
            (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
            (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
           ]

    a, b, c, d = abcd[qtype-1]
    n = len(x)
    g, j = modf( a + (n+b) * q -1)
    if j < 0:
        return y[0]
    elif j >= n:
        return y[n-1]   # oct. 8, 2010 y[n]???!! uncaught  off by 1 error!!!

    j = int(floor(j))
    if g ==  0:
       return y[j]
    else:
       return y[j] + (y[j+1]- y[j])* (c + d * g)

def calculateInterquartileRange(x):
    firstQ = quantile(x, 0.25)
    thirdQ = quantile(x, 0.75)

    return thirdQ - firstQ

def isNumeric(x):
    return str(x).replace(".", "").replace("-", "").isdigit()

def getUniqueMatrixColumnValues(filePath, columnIndex):
    uniqueValues = set()

    for line in file(filePath):
        uniqueValues.add(line.rstrip().split("\t")[columnIndex])

    return sorted(list(uniqueValues))

def fisherExactTest(x):
    return FishersExactTest.probability_of_table(x)

def complementGenomicSequence(sequence):
    mod = ""

    for base in sequence:
        mod += complementGenomicBase(base)

    return mod

def complementGenomicBase(base):
    base = base.upper()

    if base == "A":
        return "T"
    if base == "T":
        return "A"
    if base == "C":
        return "G"
    return "C"

def reverseComplementGenomicSequence(dnaSequence):
    return reverseString(complementGenomicSequence(dnaSequence))

def reverseString(string):
    return string[::-1]

def getDictValue(dictionary, key, default=""):
    if key in dictionary:
        return dictionary[key]
    return default

def getDiffPositions(string1, string2):
    matcher = difflib.SequenceMatcher(a=string1, b=string2)
    blocks = matcher.get_matching_blocks()

    diffPositions = []
    for block in blocks:
        if block[2] == 0 or block[2] == len(string1):
            continue

        if len(diffPositions) == 0:
            diffPositions.append(block[2])
        else:
            diffPositions.append(block[2] + diffPositions[-1])

    return diffPositions

def getSimilarityPercent(string1, string2):
    blocks = difflib.SequenceMatcher(None, a=string1, b=string2).get_matching_blocks()

    totalMatching = 0.0
    for block in blocks:
        totalMatching += block[2]

    return (totalMatching / float(len(string1))) * 100.0

def getLineItems(line, separator="\t"):
    return line.rstrip().split(separator)

def sortMatrix(data, columnIndex, reverse=False):
    data.sort(key=itemgetter(columnIndex), reverse=reverse)
    return data

def uniqueSort(values):
    # Slow but keeps values in order and uniquifies
    out = []

    for value in values:
        if value not in out:
            out.append(value)

    return out
