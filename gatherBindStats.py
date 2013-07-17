# Basic analysis of the diverity of the binding proteins.

import os
import sys

def getUniqueShilpaZFs(path, contacts):
	# Returns a list of ZFs for the organism of the specified
	# file (when using only given list of contact positions)
	# Assumes are in Shilpa's format and have domains listed 
	# from -1,1,2,3,4,5,6,7 positions.

	inFile = open(path, 'r')
	lines = [l.strip().split() for l in inFile if l[0] != '#']
	regions = [l[11] for l in lines if l[12] == "NOGAP"]
	prots = set()
	for r in regions:
		prot = ""
		for i in contacts:
			prot += r[i]
		prots.add(prot)

	return prots

def getProtDict(path, contacts, cut = 0):
	# Returns a dictionary mapping unique proteins 
	# (when using only given list of contact positions)
	# to a list of tuples (3-mer target bound, count).
	# path should be to a file the 'all.txt' format
	# and contatcs positions in file are -1,1,2,3,5,6
	# Can optionally remove proteins for which all
	# binding counts are less than the cut parameter.
	
	inFile = open(path, 'r')
	inFile.readline()  # Skip the header line
	protDict = {}

	for line in inFile:
		sp_line = line.strip().split('\t')
		freq = eval(sp_line[2])
		targ = sp_line[0]
		prot = ""
		for c in contacts:
			prot += sp_line[1][c]
		if freq >= cut:
			if protDict.has_key(prot):
				if protDict[prot].has_key(targ):
					protDict[prot][targ] += freq
				else:
					protDict[prot][targ] = freq
			else:
				protDict[prot] = {}
				protDict[prot][targ] = freq

	inFile.close()
	return protDict

def getProtSet(path, contacts, jsdCut):
	# Returns the set of unique proteins for the 
	# contact positions given by contacts.
	# Path should be to a file of the 'all.txt' format.
	
	inFile = open(path, 'r')
	inFile.readline()  # Skip the header line
	protSet = set()

	for line in inFile:
		sp_line = line.strip().split('\t')
		if eval(sp_line[5]) < jsdCut:
			prot = ""
			for c in contacts:
				prot += sp_line[1][c]
			protSet.add(prot)

	inFile.close()
	return protSet

def main():

	path = "../data/shilpa/Drosophila_melanogaster_ZF.fulldom"
	flyZFs = getUniqueShilpaZFs(path, [0, 1, 2, 5])
	print len(flyZFs)

if __name__ == '__main__':
	main()