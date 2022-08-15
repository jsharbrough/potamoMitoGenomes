import sys
import gzip
def seFastq2Fasta(samfile,fastq):
    readList = []
    orientDict = {}
    infile = open(samfile,'r')
    infile = open(samfile,'r')
    sys.stderr.write('Loading samfile\n')
    for line in infile:
        if line[0] != '@':
            lineSplit = line.split('\t')
            scaff = lineSplit[2]
            pos = int(lineSplit[3])
            read = lineSplit[0]
            orient = lineSplit[1]
            if read not in readList:
            	readList.append(read)
            	orientDict[read] = orient
    infile.close()
    sys.stderr.write('Samfile loaded\n')
    reads = buildFastqDict(fastq)
    sys.stderr.write('Writing to stdout\n')
    for read in readList:
        if orientDict[read] == '16' or orientDict[read] == '272' or orientDict[read] == '2064':
            currRead = reverseComplement(reads[read])
        else:    
            currRead = reads[read]
        sys.stdout.write('>' + read + '\n' + currRead + '\n')

def miseq2Fasta(samfile,refFasta,fastq1,fastq2):
    readList = []
    orientDict = {}
    refDict,refList = buildSeqDict(refFasta)
    infile = open(samfile,'r')
    fivePrimePos = len(refDict[refList[0]]) - 50
    threePrimePos = 50
    infile = open(samfile,'r')
    for line in infile:
        if line[0] != '@':
            lineSplit = line.split('\t')
            scaff = lineSplit[2]
            pos = int(lineSplit[3])
            read = lineSplit[0]
            orient = lineSplit[1]
            if scaff == refList[0] and pos >= fivePrimePos and read not in readList:
                readList.append(read)
                orientDict[read] = orient
            elif scaff == refList[1] and pos < threePrimePos and read not in readList:
                readList.append(read)
                orientDict[read] = orient
    infile.close()
    fReads = buildFastqDict(fastq1,readList)
    rReads = buildFastqDict(fastq2,readList)
    for read in readList:
        if orientDict[read] == '16' or orientDict[read] == '272' or orientDict[read] == '2064':
            fReads[read] = reverseComplement(fReads[read])
        else:    
            rReads[read] = reverseComplement(rReads[read])
        sys.stdout.write('>' + read + '_1\n' + fReads[read] + '\n>' + fivePrimePos*'-' + read + '_2\n' + fivePrimePos*'-' + rReads[read] + '\n\n')

def reverseComplement(seq):
    #seqDict = buildSeqDict(fasta)
    #for sequence in seqDict:
        #seq = seqDict[sequence]
    seq_revc = ''
    for nuc in seq:
        if nuc == 'A':
            seq_revc = 'T' + seq_revc
        elif nuc == 'T':
            seq_revc = 'A' + seq_revc
        elif nuc == 'C':
            seq_revc = 'G' + seq_revc
        elif nuc == 'G':
            seq_revc = 'C' + seq_revc
        elif nuc == 'M':
            seq_revc = 'K' + seq_revc
        elif nuc == 'R':
            seq_revc = 'Y' + seq_revc
        elif nuc == 'S':
            seq_revc = 'S' + seq_revc
        elif nuc == 'W':
            seq_revc = 'W' + seq_revc
        elif nuc == 'K':
            seq_revc = 'M' + seq_revc
        elif nuc == 'Y':
            seq_revc = 'R' + seq_revc
        elif nuc == 'V':
            seq_revc = 'B' + seq_revc
        elif nuc == 'H':
            seq_revc = 'D' + seq_revc
        elif nuc == 'D':
            seq_revc = 'H' + seq_revc
        elif nuc == 'B':
            seq_revc = 'V' + seq_revc
        else:
            seq_revc = nuc + seq_revc
    return seq_revc    
    

def buildSeqDict(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line[1:]
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict,scaffoldList

def cigarCalc(cigar):
    totalLength = 0
    currString = ''
    for char in cigar:
        if char == 'M':
            totalLength += int(currString)
            currString = ''
        elif char == 'N':
            totalLength += int(currString)
            currString = ''
        elif char == 'I':
            totalLength += int(currString)
            currString = ''
        elif char == 'D':
            currString = ''
        elif char == 'S':
            totalLength += int(currString)
            currString = ''
        elif char == 'H':
            currString = ''
        elif char == 'P':
            currString = ''
        elif char == 'X':
            totalLength += int(currString)
            currString = ''
        elif char == '=':
            totalLength += int(currString)
            currString = ''
        else:
            currString += char
    sys.stdout.write(str(totalLength))
    
def buildFastqDict(fastq):
    sys.stderr.write('Loading fq\n')
    infile = gzip.open(fastq,'r')
    seqDict = {}
    seqList = []
    lineNum = 0
    for line in infile:
        if lineNum == 0:
            seqName = line[1:].split(' ')
            seqList.append(seqName[0])
            lineNum += 1
        elif lineNum == 1:
            currSeq = line
            while currSeq[-1] == '\r' or currSeq[-1] == '\t' or currSeq[-1] == '\n':
            	currSeq = currSeq[0:-1]
            seqDict[seqName[0]] = currSeq
            lineNum += 1
        elif lineNum == 2:
            lineNum += 1
        elif lineNum == 3:
            lineNum = 0    
    sys.stderr.write('Fq Loaded\n')
    return seqDict

seFastq2Fasta(sys.argv[1],sys.argv[2])
#miseq2Fasta(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]) # miseq2Fasta(samfile,refFasta,fastq1,fastq2)