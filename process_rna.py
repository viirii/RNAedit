#!/usr/bin/env python

import sys

# How to run the script: python process_rna.py aligned_rrna.txt
inputFile = sys.argv[1]
#2d array of outputSequences and binarySequences (1 if psuedoU, else 0)
outputSeqs = []
binarySeqs = []
with open(inputFile, 'r') as iFile:
	iFileArr = iFile.readlines()
	for i in xrange(len(iFileArr)):
		#actual sequence data
		if i % 2 == 1:
			seqData = iFileArr[i]
			if "P" not in seqData:
				continue
			else:
				print i
				outputSeqStr = ""
				binarySeqStr = ""
				for n in seqData:
					if n == "-":
						continue
					elif n == "P":
						outputSeqStr += "Y"
						binarySeqStr += "1"
					else:
						outputSeqStr += n
						binarySeqStr += "0"
				outputSeqs.append(list(outputSeqStr))
				binarySeqs.append(list(binarySeqStr))

def getOutputSequences():
	return outputSeqs

def getBinarySequences():
	return binarySeqs


