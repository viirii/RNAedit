#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import sys

# How to run the script: python process_rna.py aligned_rrna.txt
inputFile = sys.argv[1]
#2d array of outputSequences and binarySequences (1 if psuedoU, else 0)
outputSeqs = []
binarySeqs = []
labels = []
symbolDict = {":": "A", "J": "U", "#": "G", "B": "C", ";": "G", "7":"G", "H":"A", 
"N": "U", "<": "C", "L": "G", "?": "C", "R": "G", "Z": "U", "M": "C", "D": "U", "=": "A",
"λ": "C", "ζ": "A",  "δ": "U", "α": "U"}
nucs = ["a", "g", "u", "c", "A", "G", "U", "C"]

def checkModNuc(n):
	if n in symbolDict: return True
	return False

with open(inputFile, 'r') as iFile:
	iFileArr = iFile.readlines()
	for i in xrange(len(iFileArr)):
		#actual sequence data
		if i % 2 == 1:
			seqData = iFileArr[i]
			if "P" not in seqData:
				continue
			else:
				outputSeqStr = ""
				binarySeqStr = ""
				for n in seqData:
					if n == "-":
						continue
					elif n == "P":
						outputSeqStr += "Y"
						binarySeqStr += "1"
					elif checkModNuc(n):
						#print n, symbolDict[n]
						outputSeqStr += symbolDict[n]
						binarySeqStr += "0"
					else:
						outputSeqStr += n
						binarySeqStr += "0"
				outputSeqs.append(outputSeqStr)
				binarySeqs.append(binarySeqStr)
				labels.append(iFileArr[i-1])

#print "outputSeqs", outputSeqs
#print "binarySeqs", binarySeqs
#print "labels", labels

def getOutputSequences():
	return outputSeqs

def getBinarySequences():
	return binarySeqs


