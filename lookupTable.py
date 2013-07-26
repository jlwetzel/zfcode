import os
import re
from pwm import *
from gatherBindStats import getProtDict

def getPosIndex(npos, canonical):
	# Return the set of indices needed depending
	# whether we need canonical sequences or 
	# complete sequences and the number of 
	# varied positions

	if npos == 6:
		if canonical:
			ind = [0, 2, 3, 5]
		else:
			ind = range(6)

	elif npos == 5:
		if canonical:
			ind = [0,1,2,4]
		else:
			ind = range(5)
	return ind

def updateTargList(fname, targList, protein, 
                   canonical, canInd):
	# Appends a tuple to targList if the protein is 
	# found in the file at path fname.
	# Assumes files named like AAA.txt	

	targ = fname.split('/')[-1].split('.')[0]
	fin = open(fname, 'r')
	totFreq = 0.0
	found = False
	# For each line of the binding file
	for line in fin:
		sp_line = line.strip().split()
		binder = sp_line[0]
		freq = eval(sp_line[1])

		# Append the single frequency for non-canonical case
		if (not canonical) and (binder == protein):
			totFreq += freq
			targList.append([targ, totFreq])
			break
		elif canonical:
			canBinder = ''
			for i in canInd:
				canBinder += binder[i]
			if canBinder == protein:
				totFreq += freq
				if not found:
					found = True

	# Append the combined frequeny for canonical case
	if canonical and found:
		targList.append([targ, totFreq])


def get3merList(dirpath, varpos, protein, canonical = False):
	# Returns a list of tuples pairs (3mer, freq),
	# where the 3 mers are DNA 3mer that bound the 
	# protein and freq is the relative frequency 
	# with which it bound.
	#
	# - dir is a path to a directory with binding files.
	# - varpos is the number of positions varied in the 
	# ZF domains in the binding assay for dir.
	# - protein is the protein (ZF) domain of interest
	# assumed to be given using the same positions varied 
	# in the binding assay by default.
	# - If canonical is set to True, then the length of 
	# protein should be 4, and domains in the binding 
	# assay will be converted canonical positions 
	# only, with the proper adjustment to freq applied.

	canInd = getPosIndex(varpos, canonical)
	targList = []  # The list of (3mer, freq) tuples to 
				   # be returned.

	# For each binding file in the dirpath
	handle = os.popen('ls ' + dirpath, 'r')
	for fname in handle:
		fname = fname.strip()
		if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
			continue

		updateTargList(dirpath + fname, targList, 
		               protein, canonical, canInd)

	# Normalize the frequencies across the bound 3mers to 1
	totFreq = 0.0
	for [targ, freq] in targList:
		totFreq += freq
	for i in targList:
		i[1] = i[1]/float(totFreq)

	return [(i[0], i[1]) for i in targList]

def targListToLogo(dstDir, targList, protein):
	# Converts a list of (3mer, frequency) tuples
	# into a sequence logo

	nucs = ['A', 'C', 'G', 'T']

	# Convert the targList into a dictionary that can
	# be passed into the 
	posCounts = initPosCounts(3, 'dna') # (pos, nuc) -> freq
	for k in targList:
		targ, freq = k[0], k[1]
		for i in range(3):
			posCounts[i, targ[i]] += freq

	pwmfile = dstDir + protein + '.txt'
	logofile = dstDir + protein + '.pdf'
	writePWM(pwmfile, posCounts, 3, nucs)

	print "Creating %s" %logofile
	makeLogo(pwmfile, logofile, alpha = 'dna', xlab = protein)

def main():
	fing = 'F2'
	strin = 'high'
	protDir = 'protein_cut10_entr04'
	npos = 6
	canonical = True
	logoDir = '/Users/jlwetzel/Desktop/F2_high_entr04_logos_lookupTable/'
	ind = getPosIndex(npos, canonical)

	dirpath = '../data/b1hData/newDatabase/6varpos/' + \
		'/'.join([fing, strin, protDir]) + '/'
	protDict = getProtDict(dirpath + 'all.txt', ind)

	bindOnly1file = open(logoDir + 'only1Target.txt', 'w')
	i = 1
	for prot in sorted(protDict.keys()):
		print "%d of %d" %(i, len(protDict))
		if len(protDict[prot]) > 1:
			targList = get3merList(dirpath, 6, prot, canonical)
			targListToLogo(logoDir, targList, prot)
		else:
			targ = protDict[prot].keys()[0]
			freq = protDict[prot][targ]
			bindOnly1file.write('%s\t%s\t%f\n' %(targ, prot, freq))
		i += 1

if __name__ == '__main__':
	main()