# Basic analysis of the diverity of the binding proteins.

import os
import sys

nucs = ['A', 'C', 'G', 'T']

def getProtDict(path, contacts):
	# Returns a dictionary mapping unique proteins 
	# (when using only given list of contact positions)
	# to a list of tuples (3-mer target bound, count).
	# path should be to a file the 'all.txt' format
	
	inFile = open(path, 'r')
	inFile.readline()  # Skip the header line
	protDict = {}

	for line in inFile:
		sp_line = line.strip().split('\t')
		count = sp_line[2]
		targ = sp_line[0]
		prot = ""
		for c in contacts:
			prot += sp_line[1][c]
		if protDict.has_key(prot):
			protDict[prot].append((targ, count))
		else:
			protDict[prot] = [(targ, count)]

	inFile.close()
	return protDict

def main():
	path = "../data/b1hData/oldData/F3/threshold025/low/all.txt"
	prots = getProtDict(path, [0, 2, 3, 5])
	print len(prots)

if __name__ == '__main__':
	main()