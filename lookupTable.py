import os
import re
import numpy as np
import math
from pwm import makeLogo, pwmfile2matrix, comparePWMs, makeNucMatFile
from fixTables import normalizeFreq
#from gatherBindStats import getProtDict

def getSubDict(fname):
	# Return a substitution dictionary indicated
	# by the given file path.

	subDict = {}

	fin = open(fname, 'r')
	xlabs = fin.readline().strip().split()
	for line in fin:
		ylab = line.strip().split()[0]
		scores = [eval(i) for i in line.strip().split()[1:]]
		for i, s in enumerate(scores):
			subDict[ylab, xlabs[i]] = s

	return subDict

###  Possible substitution matrices for nearest neighbors lookups
# WEIGHTS = None
# Use a PAM 30 matrix for weighting neighbor sequences
NEIGHBOR_WEIGHTS = getSubDict("../data/substitution_mats/PAM30.txt")

nucs = ['A', 'C', 'G', 'T']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
	      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

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

def getNeighborWeights(prot, neighbors):
	# Applies the weights to neighbors of prot 
	# based on the substituiton matrix defined 
	# by the NEIGHBOR_WEIGHTS parameter

	# Get the list of neighbor weights
	nWeights = []
	for n in neighbors:
		for i, a in enumerate(n):
			if prot[i] != a:
				nScore = math.exp(NEIGHBOR_WEIGHTS[prot[i], a])
		nWeights.append(nScore)

	# Shift so that min weight is zero then normalize 
	# all weights to between 0 and 1
	nWeights = np.array(nWeights, dtype=float)
	nWeights = nWeights - np.min(nWeights)
	nWeights = nWeights/np.sum(nWeights)
	return nWeights

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

def updateTargListNN(fname, targList, protein, 
                     canonical, canInd, neighbors,
                     nWeights = None):
	# Appends a tuple to targList if a nearest
	# neighbor to the protein is found when 
	# allowing the 1, 3, or 6 positions of 
	# the alpha helix to vary.
	# Assumes files named like AAA.txt	

	targ = fname.split('/')[-1].split('.')[0]
	fin = open(fname, 'r')
	found = False
	
	# For each line of the binding file
	totFreq = 0.0
	for line in fin:
		sp_line = line.strip().split()
		freq = eval(sp_line[1])

		# Get the binding protein for this line
		if not canonical:
			binder = sp_line[0]
		else:
			binder = ''
			for i in canInd:
				binder += sp_line[0][i]
		
		# Add to the frequency if this binder is a neighbor,
		# correcting for weighting if necessary
		for i, n in enumerate(neighbors):
			if binder == n:
				if NEIGHBOR_WEIGHTS != None:
					totFreq += freq*nWeights[i]
				else:
					totFreq += freq
				if not found:
					found = True
					break

	# Append to the targList
	if found:
		targList.append([targ, totFreq])

def decomposeNeighbors(protein, neighbors, decompose):
	# Returns a dictionary of neighbor lists according 
	# to the decompose dictionary
	#
	# protein is a protein sequence
	# neighbors is the list of "off-by-one" neighbors for protein
	# decompose is a dictionary, indexed by base positions,
	# where each index points to a list of protein positions not 
	# allowed to vary for the neighbor list for that index

	neighborDict = {}
	# For each base position
	for bpos in decompose.keys():
		neighborDict[bpos] = []
		# For each neighbor
		for n in neighbors:

			# For each amino position to remain fixed for this bpos
			keepNeighbor = True
			for apos in decompose[bpos]:
				keepNeighbor = keepNeighbor and (protein[apos] == n[apos])
			if keepNeighbor:
				neighborDict[bpos].append(n)

	return neighborDict

def normalizeTargList(targList):
	# targList is a list of tuples where first element 
	# is a target and second element is a frequency with 
	# which the target was bound.  This function just 
	# normalizes the distribution of frequencies.
	totFreq = 0.0
	for [targ, freq] in targList:
		totFreq += freq
	for i in targList:
		i[1] = i[1]/float(totFreq)

def get3merList(dirpath, varpos, protein, canonical = False,
				useNN = False, skipExact = False, decompose = None):
	# Returns a list of tuples pairs (3mer, freq),
	# where the 3 mers are DNA 3mer that bound the 
	# protein and freq is the relative frequency 
	# with which it bound.
	#
	# - dir is a path to a directory with binding files.
	#
	# - varpos is the number of positions varied in the 
	# ZF domains in the binding assay for dir.
	#
	# - protein is the protein (ZF) domain of interest
	# assumed to be given using the same positions varied 
	# in the binding assay by default.
	#
	# - If canonical is set to True, then the length of 
	# protein should be 4, and domains in the binding 
	# assay will be converted canonical positions 
	# only, with the proper adjustment to freq applied.
	#
	# - If useNN, then look for nearest neighbors when can't
	# find an exact match in the binding assay
	# (see updateTargListNN for details)
	#
	# - If skipExact, then skip the search for exact matches
	# and go straight to nearest neighbors

	canInd = getPosIndex(varpos, canonical)
	targList = []  # The list of (3mer, freq) tuples to 
				   # be returned.

	# For each binding file in the dirpath, look for exact 
	# matches first.
	handle = os.popen('ls ' + dirpath, 'r')
	for fname in handle:
		if skipExact:
			break
		fname = fname.strip()
		if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
			continue

		updateTargList(dirpath + fname, targList, 
		               protein, canonical, canInd)

	# If the list is empty after looking for exact matches,
	# then look for nearest neighbors if useNN flag is set.
	if useNN:
		# Get the list of possible neighboring seqs
		neighbors = []
		if canonical:
			#swapInd = [0,2,3]  # old way when fixing position 2
			swapInd = range(4)
		else:
			#swapInd = [0,3,5]  # old way when fixing pos. 2
			swapInd = [1, 2, 3, 5]		
		for i in swapInd:
			for a in aminos:
				if a != protein[i]:
					neighbors.append(protein[:i]+a+\
					                 protein[i+1:])

		if decompose == None:
			# Get weights for the potential neighbors based on 
			# a substitution matrix
			if NEIGHBOR_WEIGHTS != None:
				nWeights = getNeighborWeights(protein, neighbors)

			# Get the targetList for this protein
			handle = os.popen('ls ' + dirpath, 'r')
			for fname in handle:
				fname = fname.strip()
				if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
					continue
				updateTargListNN(dirpath + fname, targList, 
		        	             protein, canonical, canInd, 
		            	         neighbors, nWeights)

			# Normalize the frequencies across the bound 3mers to 1
			normalizeTargList(targList)
			return targList
			
		else:
			# Decompose the list of neghbors on a per base position 
			# basis, dictated by the dictionary "decompose", and 
			# return a dictionary of targLists, indexed by base position.
			neighborDict = decomposeNeighbors(protein, neighbors, decompose)
			#print decompose
			#for k in neighborDict.keys():
			#	print "%s: %s -- %d" %(k, repr(neighborDict[k]), len(neighborDict[k]))
			targListDict = {}
			for k in neighborDict.keys():
				targListDict[k] = []

				# Get weights for each neighbor
				if NEIGHBOR_WEIGHTS != None:
					nWeights = getNeighborWeights(protein, neighbors)
				
				# Get the target list for this base
				handle = os.popen('ls ' + dirpath, 'r')
				for fname in handle:
					fname = fname.strip()
					if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
						continue
					updateTargListNN(dirpath + fname, targListDict[k], 
			        	             protein, canonical, canInd, 
			            	         neighborDict[k], nWeights)
				normalizeTargList(targListDict[k])

			#for k in targListDict:
			#	print "%s: %s -- %d" %(k, repr([i for i in targListDict[k]]), 
			#	                       len(targListDict[k]))

			return targListDict


def targListToFreqMat(targList):
	# Converts the 3-mer list into of form
	# where each element is a tuple of the form
	# (3mer, freq with which protein is seen)
	# to a 2D position frequency matrix.
	# Returns the frequency matrix

	pwm = np.zeros((3,4), float)
	for i in range(len(targList)):
		targ = targList[i][0]
		freq = targList[i][1]
		for j, n in enumerate(targ):
			pwm[j,nucs.index(n)] += freq

	return pwm

def makeDir(path):
	# Try to make a path and pass if it can't be created
	try:
		os.mkdir(path)
	except OSError:
		pass

def lookupCanonZFArray(inDir, canonZFs, useNN = True, 
                       skipExact = False, decompose = None):
	# Performs modular lookup for an array of canonical 
	# helix-position ZF domains.  Domains should be given
	# in reverse order of their appearance on the protein.

	numZFs = len(canonZFs)
	pwm = np.zeros(shape = (numZFs*3, 4))
	
	for i in range(numZFs):
		nmat = lookupCanonZF(inDir, canonZFs[i], 
		                     useNN, skipExact,
		                     decompose)
		for j in range(len(nmat)):
			pwm[i*len(nmat) + j,:] = nmat[j,:]
	
	return pwm

def singleColTargListToFreqVector(bpos, targList):
	# Returns the frequency vector corresponding to the 
	# the frequency of A, C, G, T in for base bpos
	#
	# targList is a list of two elecment lists, where the 
	# first element of each inner list should be a nuc 3mer
	# and the second element should be the frequency with 
	# which a protein of interest was bound by the 3mer.

	bposind = bpos - 1
	baseFreqList = [0.0, 0.0, 0.0, 0.0]
	for i in targList:
		base = i[0][bposind]
		freq = i[1]
		baseFreqList[nucs.index(base)] += freq

	vect = np.array(baseFreqList, dtype=float)
	vect = vect/np.sum(vect)
	#print targList
	#print vect
	return vect

def lookupCanonZF(inDir, canonProt, useNN = True, skipExact = False,
                  decompose = None):
	# Takes as input a ZF domain (canoical positions  only,
	# helix positions -1, 2, 3, 6) and returns the predicted 
	# 3-base binding specificity as a normalized 2d frequency 
	# matrix.  The matrix is produced by using the 
	# reverse look-up table approach on the forward selection 
	# data directory inDir.
	#
	# Is useNN is set to True, then then if we can't find the 
	# domain by direct reverse lookup, we look at nearest neighbors
	# by varying each of helix positions (-1,3,6) one at a time.
	# If skipExact is True, then we skip the direct look-up and 
	# go directly to nearest neighbors.
	#
	# If decompose is not set to None, it should be set to a 
	# dictionary mapping each base position to a list of important
	# amino positions for predicting that base.  These preferences
	# will then be used to decompose the nearest neighbor lookup.


	# B1H forward experiments.  Need to update if start 
	# using the 5 position experiments.
	npos = 6
	canonical = True
	ind = getPosIndex(npos, canonical)

	if decompose == None:
		# Make the list of targets bound and normalized frequencies.
		targList = get3merList(inDir, 6, canonProt, canonical,
		                       useNN, skipExact, decompose)

		# If the target list is empty, we can't return a specificity,
		# so we just return a uniform distribution for each base
		if targList == []:
			nucmat = np.zeros((3,4), float)
			for i in range(len(nucmat)):
				for j in range(len(nucmat[i,:])):
					nucmat[i,j] = 0.25
			return nucmat

		# Convert the target list to a position freq mat and return
		return targListToFreqMat(targList)

	else:
		# Get a dictionary mapping each base position to its own 
		# specific target list based on the decompose dictionary
		targLists = get3merList(inDir, 6, canonProt, canonical,
		                        useNN, skipExact, decompose)
		nucmat = np.zeros((3,4), float)
		for i, k in enumerate(sorted(targLists.keys())):
			if targLists[k] == []:
				nucmat[i,:] = np.array([0.25, 0.25, 0.25, 0.25])
			else:
				nucmat[i,:] = singleColTargListToFreqVector(k, targLists[k])
		return nucmat

def lookupMarcusPWMs(inDir, outputDir, finger, strin, 
                     filt, pred, 
                     useNN = True, skipExact = False,
                     decompose = None):
	# Make predcitions for each of the proteins that 
	# Marcus made experimental PWMs for

	# Create the prediction directory structure
	predictionDir = outputDir+'predictions/'
	pwmdir = predictionDir + 'pwms/'
	logodir = predictionDir + 'logos/'
	makeDir(predictionDir)
	makeDir(pwmdir)
	makeDir(logodir)
	
	# Get the directory containing the 3 position
	# experimental PWMs for this finger and set up
	# the output file for writing results
	if finger == 'F2':
		expDir = '../data/revExp/F2_GAG/pwms3/'
	elif finger == 'F3':
		expDir = '../data/revExp/F3_GCG/pwms3/'
	fout = open(predictionDir + 'compare.txt', 'w')
	# Write header to results file
	fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
	           %('num', 'targ','prot','canonprot','score','colcor', \
	           'colcorIC', 'totcol', 'pred.filt', 'finger', 'strin'))
	
	
	# B1H forward experiments.  Need to update if start 
	# using the 5 position ones
	npos = 6
	canonical = True
	ind = getPosIndex(npos, canonical)
	
	# Try to do a lookup table prediction for each of the 
	# proteins for which we have an experimental 3pos pwm
	for fname in os.popen('ls ' + expDir):

		# Skip if this is the wrong stringency.
		fname = fname.strip()
		if strin == 'low' and re.match(r'(.)*_5mM.txt', fname) == None:
			continue
		if strin == 'high' and re.match(r'(.)*_20mM.txt', fname) == None:
			continue

		# Get the info about this experiment
		sp_fname = fname.split('_')
		protNum = sp_fname[0]
		goal = sp_fname[1]
		prot = sp_fname[2].split('.')[0]
		canonProt = prot[0] + prot[2] + prot[3] + prot[6]
		label = '_'.join([str(protNum), goal, prot, strin])
	
		if decompose == None:
			# Make the list of targets bound and normalized frequencies.
			targList = get3merList(inDir, 6, canonProt, canonical,
		    	                   useNN, skipExact, decompose)

			# Do something else if the target list is empty
			if targList == []:
				continue

			# Convert the 3mer target list to a frequency matrix,
			# write to file, and create a logo
			nucMat = targListToFreqMat(targList)
		
		else:
			#print protNum, canonProt
			# Get a dictionary mapping each base position to its own 
			# specific target list based on the decompose dictionary
			targLists = get3merList(inDir, 6, canonProt, canonical,
			                        useNN, skipExact, decompose)
			nucMat = np.zeros((3,4), float)
			for i, k in enumerate(sorted(targLists.keys())):
				if targLists[k] == []:
					nucMat[i,:] = np.array([0.25, 0.25, 0.25, 0.25])
				else:
					nucMat[i,:] = singleColTargListToFreqVector(k, targLists[k])


		makeNucMatFile(pwmdir, label, nucMat)
		logoIn = pwmdir + label + '.txt'
		logoOut = logodir + label + '.pdf'
		makeLogo(logoIn, logoOut, alpha = 'dna', 
		         colScheme = 'classic',
		         annot = "'5,M,3'",
		         xlab = '_'.join([goal,prot]))
		
		# Compare this pwm to the reverse experiment
		score, colcor, colcorIC, totCol = \
			comparePWMs(nucMat, pwmfile2matrix(expDir + fname))

		# Write the comparison results to file
		fout.write("%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\t%s\t%s\t%s\n" \
		           %(protNum, goal, prot, canonProt, score, colcor, \
		           	 colcorIC, totCol, pred+'.'+filt, finger, strin))
	fout.close()

def main():

	"""
	# Don't use nearest neighbors
	fings = ['F2']
	strins = ['low']
	filts = ['cut10bc_0_5', 'cut3bc_0_5'] #[ 'cut10bc_0', 'cut3bc_025', 'cut10bc_025']
	filtsLabs = ['c10_0_5', 'c3_0_5'] #['c3_025', 'c10_025', 'c10_0']
	for f in fings:
		for s in strins:
			for i, filt in enumerate(filts):
				inDir = '../data/b1hData/newDatabase/6varpos/' \
					+ f + '/' + s + '/' + 'protein_seq_' + filt + '/'
				outDir = '../data/lookupTable/'	+ f + '/' + s + \
					'/' + filt + '/'
				
				lookupMarcusPWMs(inDir, outDir, f, s, filtsLabs[i],
				                 'look', useNN = False, skipExact = False)

	# Use nearest neighbors if exact matches can't be found
	fings = ['F2']
	strins = ['low']
	filts = ['cut3bc_025', 'cut10bc_025', 'cut10bc_0']
	filtsLabs = ['c3_025', 'c10_025', 'c10_0']
	for f in fings:
		for s in strins:
			for i, filt in enumerate(filts):
				inDir = '../data/b1hData/newDatabase/6varpos/' \
					+ f + '/' + s + '/' + 'protein_seq_' + filt + '/'
				outDir = '../data/lookupTableNN/'	+ f + '/' + s + \
					'/' + filt + '/'
				
				lookupMarcusPWMs(inDir, outDir, f, s, filtsLabs[i],
				                 'look', useNN = TRUE, skipExact = False)
	"""
	decomp1 = {1: [3], 2: [2,3], 3: [0,1]}

	# Use nearest neighbors and skip all exact matches
	fings = ['F2']
	strins = ['low']
	filts = ['cut10bc_0_5', 'cut3bc_0_5', 'cut10bc_0', 'cut3bc_025', 'cut10bc_025']
	filtsLabs = ['c10_0_5', 'c3_0_5', 'c3_025', 'c10_025', 'c10_0']
	for f in fings:
		for s in strins:
			for i, filt in enumerate(filts):
				#print f, s, filt
				inDir = '../data/b1hData/newDatabase/6varpos/' \
					+ f + '/' + s + '/' + 'protein_seq_' + filt + '/'
				outDir = '../data/lookupTableNNonly_pam30_decomp1/' + f + '/' + s + \
					'/' + filt + '/'
				lookupMarcusPWMs(inDir, outDir, f, s, filtsLabs[i],
				                 'look.nnOnly.PAM30_decomp1', useNN = True, skipExact = True,
				                 decompose = decomp1)


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