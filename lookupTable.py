import os
import re
from pwm import *
from gatherBindStats import getProtDict
from fixTables import normalizeFreq

nucs = ['A', 'C', 'G', 'T']

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

def targListToLogo(dstDir, targList, protein,
                   targ, label):
	# Converts a list of (3mer, frequency) tuples
	# into a sequence logo

	# Convert the targList into a dictionary that can
	# be passed into the 
	posCounts = initPosCounts(3, 'dna') # (pos, nuc) -> freq
	for k in targList:
		targ, freq = k[0], k[1]
		for i in range(3):
			posCounts[i, targ[i]] += freq
	###
	#if protein == 'LNDHLQN':
	#	print posCounts

	pwmfile = dstDir + 'pwms/' + label + '.txt'
	logofile = dstDir + 'logos/' + label + '.pdf'
	writePWM(pwmfile, posCounts, 3, nucs)

	print "Creating %s" %logofile
	makeLogo(pwmfile, logofile, alpha = 'dna', 
	         colScheme = 'classic',
	         annot = "'5,M,3'",
	         xlab = '_'.join([targ,protein]))

def lookup700s(inDir, outputDir):
	# Make predictions for the F2 reverse experiments.

	predictionDir = outputDir+'predictions/'
	expDir = '../data/revExp/F2_GAG/pwms3/'
	fout = open(predictionDir + 'F2compare.txt', 'w')
	fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %('num', \
	           'targ','prot','canonprot','score','colcor', \
	           'colcorIC', 'totcol'))
	
	#npos = 6
	#canonical = True
	#ind = getPosIndex(npos, canonical)
	#protDict = getProtDict(inDir + 'all.txt', ind)

	for fname in os.popen('ls ' + expDir):
		# Get the info about this prediction
		fname = fname.strip()
		sp_fname = fname.split('_')
		protNum = sp_fname[0]
		targ = sp_fname[1]
		prot = sp_fname[2].split('.')[0]
		canonProt = prot[0] + prot[2] + prot[3] + prot[6]
		label = '_'.join([str(protNum), targ, prot])
		
		# Make the prediciton and write out to the 
		# correct files.
		targList = get3merList(inDir, 6, canonProt, canonical)
		if targList == []:
			continue
			# Apply nearest neighbor strategy here if targList 
			# is empty instead of just ignoring?
		targListToLogo(predictionDir, targList, prot, targ, label)
		
		# Compare this pwm to the reverse experiment
		score, colcor, colcorIC, totCol = \
			comparePWMs(pwmfile2matrix(predictionDir + \
			            'pwms/' + label + '.txt'), 
		                pwmfile2matrix(expDir + label + '.txt'))
		fout.write("%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\n" %(protNum, \
		           targ, prot, canonProt, score, colcor, \
		           colcorIC, totCol))
	fout.close()

def main():


	inDir = '../data/b1hData/newDatabase/6varpos/F2/low/' + \
		'protein_seq_cut10bc_0/'
	outDir = '../data/lookupTable/cut10bc_0/'
	lookup700s(inDir, outDir)


if __name__ == '__main__':
	main()


## Old version fo main that predicts for all extant sequences
## in a given protein directory

"""
def main():
	fing = 'F2'
	strin = 'high'
	protDir = 'threshold025'
	npos = 6
	canonical = True
	logoDir = '/Users/jlwetzel/Desktop/Anton_F2_high_entr025' + \
		'_logos_lookupTable/'
	ind = getPosIndex(npos, canonical)
	dirStyle = 'old'    # For Anton's files or my files


	if dirStyle == 'new':
		dirpath = '../data/b1hData/newDatabase/6varpos/' + \
			'/'.join([fing, strin, protDir]) + '/'
	elif dirStyle == 'old':
		dirpath = '../data/b1hData/oldData/' + \
			'/'.join([fing, protDir, strin]) + '/'

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
			bindOnly1file.write('%s\t%s\t%f\n' %(prot, targ, freq))
		i += 1
		"""