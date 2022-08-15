import sys
def cigarCalc(cigar):
	totalLength = 0
	softClipped = 0
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
	return totalLength


def mergeReads(read):
	pairedRead = readDict[read]
	r1 = pairedRead[0]
	r2 = pairedRead[1]
	seq1 = r1['seq']
	seq2 = r2['seq']
	end = r1['end']
	pos2 = r2['pos2']
	if pos2 <= end:
		ov = end-pos2
		sr = seq1[0:-ov] + seq2
	else:
		sr = seq1 + 'N'*(pos2-end) + seq2
	sys.stdout.write('>' + read + '_super-read\n' + sr + '\n')

def trimSoftClip(cigar):
	softClipped = 0
	currString = ''
	for char in cigar:
		if char == 'M':
			break
			currString = ''
		elif char == 'N':
			currString = ''
		elif char == 'I':
			currString = ''
		elif char == 'D':
			currString = ''
		elif char == 'S':
			softClipped += int(currString)
			currString = ''
		elif char == 'H':
			currString = ''
		elif char == 'P':
			currString = ''
		elif char == 'X':
			break
			currString = ''
		elif char == '=':
			break
			currString = ''
		else:
			currString += char
	return softClipped
readDict = {}
readList = []
samfile = open(sys.argv[1],'r')
for line in samfile:
	realLine = line
	while realLine[-1] == '\t' or realLine[-1] == '\r' or realLine[-1] == '\n':
		realLine = realLine[0:-1]
	lineSplit = realLine.split('\t')
	pos = int(lineSplit[3])
	length = cigarCalc(lineSplit[5])
	read = lineSplit[9]
	if lineSplit[0] not in readDict:
		end = pos + length
		r1 = {'seq':read,'end':end}
		readDict[lineSplit[0]] = [r1]
		readList.append(lineSplit[0])
	elif len(readDict[lineSplit[0]]) == 1:
		r2 = {'seq':read[trimSoftClip(lineSplit[5]):],'pos2':pos}
		readDict[lineSplit[0]] += [r2]
	
for read in readList:
	mergeReads(read)