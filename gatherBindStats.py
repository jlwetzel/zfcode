# Basic analysis of the diverity of the binding proteins.

import os
import sys
import re

nucs = ['A', 'C', 'G', 'T']

def getUniqueShilpaZFs(path, contacts):
	# Returns a list of ZFs for the organism of the specified
	# file (when using only given list of contact positions)
	# Assumes are in Shilpa's format and have domains listed 
	# from -1,1,2,3,4,5,6,7 positions.

	inFile = open(path, 'r')
	#C2H2regex = re.compile( \
	#    'C[A-Z]{2,5}C[A-Z]{3,3}[FY][A-Z]{7,7}H[A-Z]{3,5}[HC]')
	lines = [l.strip().split() for l in inFile if l[0] != '#']
	regions = [l[-1] for l in lines \
		if l[12] == "NOGAP" and [6] != "zf-C2H2_4"]# and \
		#C2H2regex.search(l[-1]) != None]
	prots = set()
	for r in regions:
		prot = ""
		for i in contacts:
			prot += r[i]
		prots.add(prot)

	return prots

def getProtDict(path, contacts):
	# Returns a dictionary mapping unique proteins 
	# (when using only given list of contact positions)
	# to a list of tuples (3-mer target bound, count).
	# path should be to a file the 'all.txt' format
	# and contatcs positions in file are -1,1,2,3,5,6
	
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
	path = "../data/shilpa/Drosophila_melanogaster_ZF.fulldom"
	flyZFs = getUniqueShilpaZFs(path, [0, 1, 2, 5])
	print len(flyZFs)

if __name__ == '__main__':
	main()