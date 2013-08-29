import os
import re
import numpy as np
import math
from pwm import makeLogo, pwmfile2matrix, comparePCC, getConsensus, makeNucMatFile
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
NEIGHBOR_WEIGHTS = getSubDict('../data/substitution_mats/PAM30.txt')
NEIGHBOR_ORDER = getSubDict('../data/substitution_mats/PAM30.txt')

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

def getNeighborWeights(prot, neighbors, skipExact, matType = 'weight'):
	# Returns a numpy array of weights for each neighbor in the 
	# order that the neighbors appear in the list.
	# matType decides which matrix to use, either the 
	# one for ordering or the one for weighting.

	# Decide which matrix to use
	if matType == 'weight':
		weightMat = NEIGHBOR_WEIGHTS
	elif matType == 'order':
		weightMat = NEIGHBOR_ORDER

	exactMatchInList = False

	# Get the list of neighbor weights
	nWeights = []
	for k, n in enumerate(neighbors):

		# Save a spot for the exact match if it's in this list
		if not skipExact and n == prot:
			exactMatchInList = True
			exactMatchInd = k
			nWeights.append(0.0)
		
		# Otherwise give the score of the exponential 
		# from the substitution matrix for the amino substitution
		else:
			for i, a in enumerate(n):
				if prot[i] != a:
					nScore = math.exp(weightMat[prot[i], a])
			nWeights.append(nScore)

	# If used exact matching neighbor, then make its 
	# weight as heavy as the heaviest non-exact neighbor
	nWeights = np.array(nWeights, dtype=float)
	if not skipExact and exactMatchInList:
		nWeights[exactMatchInd] = np.max(nWeights)

	# Normalize weights to a distribution
	nWeights = nWeights/np.sum(nWeights)

	return nWeights

def decomposeNeighbors(protein, neighbors, decompose, skipExact):
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

			# Include the exact match for each base position if not
			# skipping exact matches
			if not skipExact and n == protein:
				neighborDict[bpos].append(n)
				continue

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
	return vect

def getNeighborBaseVect(freqDict, neighbor, bpos):
	# Returns the vector of base frequencies for 
	# bpos for all matches to this neighboring protein

	# Get a target list for this neighbor
	targList = []
	for targ in freqDict.keys():
		if freqDict[targ].has_key(neighbor):
			targList.append([targ, freqDict[targ][neighbor]])

	# In the case that neighbor bound nothing, return None
	if targList == []:
		return None

	# Otherwise convert the targList to a normalized vector 
	# of frequencies for the base position indicated
	else:
		normalizeTargList(targList)
		return singleColTargListToFreqVector(bpos, targList)

def computeFreqDict(dirpath, ind): 
	# Run through all the files for this binding dictionary 
	# and create one large binding dictionary of the form:
	# 3mer -> bindingProtein -> freq

	freqDict = {}
	handle = os.popen('ls ' + dirpath, 'r')
	for fname in handle:
		fname = fname.strip()
		if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
			continue
		targ = fname.split('.')[0]
		freqDict[targ] = {}

		fin = open(dirpath + fname, 'r')
		for line in fin:
			sp_line = line.strip().split()
			binder = sp_line[0]
			freq = eval(sp_line[1])
			
			bindSub = ''
			for i in ind:
				bindSub += binder[i]

			if freqDict[targ].has_key(bindSub):
				freqDict[targ][bindSub] += freq
			else:
				freqDict[targ][bindSub] = freq

	return freqDict

def sortAndWeightNeighbors(prot, bpos, neighbors, skipExact):
	# Sorts the neighbors in the order in which we 
	# would like to process them (most "important"
	# neighbors first).
	# 
	# Returns the newly sorted list of neighbors and 
	# a list of weights in the same order.
	
	# The order in which aminos should be allowed to vary 
	# for each of the bases.
	order = {1: [0, 1, 2],  # -1,2,3
			 2: [1, 0, 3],  #  2,-1,6
			 3: [3, 2, 1]}  #  6,3,2

	# Separate neighbors according to varied positions
	neighborBins = {}
	for n in neighbors:
		
		# Ignore the exact match for now if using it
		if not skipExact and n == prot:
			continue

		for i, a in enumerate(n):
			if a != prot[i]:
				varpos = i
				break
		if neighborBins.has_key(varpos):
			neighborBins[varpos].append(n)
		else:
			neighborBins[varpos] = [n]

	# If using the exact neighbor, add it now to the 
	# first bin
	if not skipExact:
		neighborBins[order[bpos][0]].append(prot)
	
	# Sort within each neighbor group according 
	# to the matrix associated with NEIGHBOR_ORDER variable
	neighborsSorted = []
	for k in order[bpos]:
		if NEIGHBOR_ORDER != None:
			order_weights = getNeighborWeights(prot, neighborBins[k], 
			                                   skipExact, matType = 'order')
		else:
			unifWeight = 1/float(len(neighborBins[k]))
			order_weights = np.array([unifWeight]*len(neighborBins[k]),
			                         dtype = float)
		# Get the indices for the reverse sorted order
		orderInd = [i[0] for i in sorted(enumerate(order_weights),
		                                 key = lambda x:x[1],
		                                 reverse = True)]

		# Add this bin's neighbors to the list in the sorted order
		for i in orderInd:
			neighborsSorted.append(neighborBins[k][i])
	
	# Get the weights for the now sorted complete list of neighbors
	if NEIGHBOR_WEIGHTS != None:
		nWeights = getNeighborWeights(prot, neighborsSorted,
		                              skipExact, matType = 'weight')
	else:
		unifWeight = 1/float(len(neighborsSorted))
		nWeights = np.array([unifWeight]*len(neighborsSorted),
			                     dtype = float)

	return neighborsSorted, nWeights

def getTopKNeighborsPWM(freqDict, prot, neighborDict, topk, skipExact):
	# Returns a numpy array pwm of for the prediciton
	# given the top k neighbors only in terms of distance
	# from the original protein.  The neighbor distances 
	# are ranked in a hierarchical fashion.  First level 
	# of the hierarchy is based on the position which has 
	# been allowed to vary.  Position farthest form the
	# canonical base partner are given highest closeness.
	# Within each position, we choose in order of the 
	# rankings by the weight matrix currently being 
	# used.
	# 
	# If topk is set to None, then us all neighbors
	# ... otherwise use the topk neighbors
	#
	# neighborDict is a dictionary of neighbors indexed by 
	# base.  E.g. neighborDict[1] points to all alowable 
	# neighbors for base position 1

	# A numpy array to be returned after filling in
	pwm = np.zeros((3,4), dtype = float)
	
	# For each base
	# print neighborDict.keys()
	neighborsPerBase = []
	for k in neighborDict.keys():
		baseVectors = []
		neighborsUsed = 0
		
		# Sort the neighbors in the order in which we'd like 
		# to use them.
		if topk != None:
			#print k
			neighborDict[k], nWeights = sortAndWeightNeighbors(prot, k,
			                                                   neighborDict[k],
			                                                   skipExact)
		else:
			if NEIGHBOR_WEIGHTS != None:
				nWeights = getNeighborWeights(prot, neighborDict[k], skipExact)
			else:
				unifWeight = 1/float(len(neighborDict[k]))
				nWeights = np.array([unifWeight]*len(neighborDict[k]),
			                       dtype = float)

		# Get the normalized frequency vectors for each neighbor at  
		# this base position
		#print
		#print len(neighborDict[k])
		for i, n in enumerate(neighborDict[k]):
			newVect = getNeighborBaseVect(freqDict, n, k)
			if newVect != None:
				pass
				#print "Base: %d\tNeighbor: %s\tWeight: %.5f" %(k, n, nWeights[i])
				#print newVect
			baseVectors.append(newVect)
		
		# For each neighbor found, weight its vector by that
		# neighbor's weight.
		for i in range(len(nWeights)):
			if baseVectors[i] != None:
				baseVectors[i] = baseVectors[i] * nWeights[i]

		# Remove Nones from the baseVector list and 
		# record how many neighbors were used
		baseVectors = [i for i in baseVectors if i != None]
		

		# If we found at least one neighbor
		if baseVectors != []:
			# Add the weighted vectors together
			for i in range(len(baseVectors)):
				if neighborsUsed == topk:
					break
				else:
					pwm[k-1] = pwm[k-1] + baseVectors[i]
					neighborsUsed += 1
			# Renormalize since some neighbors may not have been used
			pwm[k-1] = pwm[k-1]/np.sum(pwm[k-1]) 
			neighborsPerBase.append(neighborsUsed)
			#print "Used %d neighbors for base %d" %(neighborsUsed, k)
		
		# If we found no neighbors use a uniform vector
		else:
			pwm[k-1] = np.array([0.25, 0.25, 0.25, 0.25], dtype = float)

	#print neighborsPerBase
	return pwm, neighborsPerBase


def get3merList(freqDict, protein, canonical = False,
				useNN = False, skipExact = False, 
				decompose = None, topk = None):
	# Returns a list of tuples pairs (3mer, freq),
	# where the 3 mers are DNA 3mer that bound the 
	# protein and freq is the relative frequency 
	# with which it bound.
	#
	# - dir is a path to a directory with binding files.
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

	# Read the entire list of target/binder 
	# frequencies into a hash table

	targList = []  # The list of (3mer, freq) tuples to 
				   # be returned.

	# If not using nearest neighbors, look for exact matches
	# in the lookup table
	if not useNN:
		for targ in sorted(freqDict.keys()):
			if skipExact:
				break
			if freqDict[targ].has_key(protein):
				targList.append([targ, freqDict[targ][protein]]) 
		normalizeTargList(targList)
		print "Here"
		return targList

	# The target list is empty and we are allowed to use 
	# nearest neighbors
	neighbors = []
	if canonical:
		swapInd = range(4)
	else:
		swapInd = [0, 2, 3, 5]		
	for i in swapInd:
		for a in aminos:
			if a != protein[i]:
				neighbors.append(protein[:i]+a+\
				                 protein[i+1:])
	# Include the exact match as a neighbor if not told otherwise
	if not skipExact:
		neighbors.append(protein)

	# Get the per-base neighbor decomposition
	neighborDict = decomposeNeighbors(protein, neighbors, decompose, skipExact)
	
	# Return the pwm obtained by decomposing neighbors
	# and using the top k of them

	# Here

	return getTopKNeighborsPWM(freqDict, protein, neighborDict, topk, skipExact)
	
def makeDir(path):
	# Try to make a path and pass if it can't be created
	try:
		os.mkdir(path)
	except OSError:
		pass

def lookupCanonZFArray(freqDict, canonZFs, useNN = True, 
                       skipExact = False, decompose = None, topk = None):
	# Performs modular lookup for an array of canonical 
	# helix-position ZF domains.  Domains should be given
	# in reverse order of their appearance on the protein.

	numZFs = len(canonZFs)
	pwm = np.zeros(shape = (numZFs*3, 4))
	
	for i in range(numZFs):
		nmat = lookupCanonZF(freqDict, canonZFs[i], 
		                     useNN, skipExact, decompose, topk)
		for j in range(len(nmat)):
			pwm[i*len(nmat) + j,:] = nmat[j,:]
	
	return pwm

def lookupCanonZF(freqDict, canonProt, useNN = True, skipExact = False,
                  decompose = None, topk = None):
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

	# Make the list of targets bound and normalized frequencies.
	if not useNN:
		targList = get3merList(freqDict, canonProt, canonical,
	                           useNN, skipExact, decompose, topk)
	else:
		targList, neighborsPerBase = get3merList(freqDict, canonProt, canonical,
	                                         useNN, skipExact, decompose, topk)
	
	# This is target/frequncy list
	if isinstance(targList, list):
		# If the target list is empty, we can't return a specificity,
		# so we just return a uniform distribution for each base
		print "Here"
		if targList == []:
			nucmat = np.zeros((3,4), float)
			for i in range(len(nucmat)):
				for j in range(len(nucmat[i,:])):
					nucmat[i,j] = 0.25
			return nucmat

		# Convert the target list to a position freq mat and return
		return targListToFreqMat(targList)

	# This is already a proper numpy array
	else:
		return targList, neighborsPerBase

def setWeightMatrices(weight_mat, order_mat):
	# Set the global weight matrix parameters
	global NEIGHBOR_WEIGHTS
	global NEIGHBOR_ORDER
	
	if weight_mat == 'PAM30':
		NEIGHBOR_WEIGHTS = getSubDict("../data/substitution_mats/PAM30.txt")
	elif weight_mat == 'PAM100':
		NEIGHBOR_WEIGHTS = getSubDict("../data/substitution_mats/PAM100.txt")
	elif weight_mat == 'PAM250':
		NEIGHBOR_WEIGHTS = getSubDict("../data/substitution_mats/PAM250.txt")
	else:
		NEIGHBOR_WEIGHTS = None

	if order_mat == 'PAM30':
		NEIGHBOR_ORDER = getSubDict("../data/substitution_mats/PAM30.txt")
	elif order_mat == 'PAM100':
		NEIGHBOR_ORDER = getSubDict("../data/substitution_mats/PAM100.txt")
	elif order_mat == 'PAM250':
		NEIGHBOR_ORDER = getSubDict("../data/substitution_mats/PAM250.txt")
	else:
		NEIGHBOR_ORDER = None

def main():

	# Debugging stuff
	inDir = '../data/b1hData/newDatabase/6varpos/F2/low/protein_seq_cut10bc_0_5/'
	canonical = True
	varpos = 6
	canInd = getPosIndex(varpos, canonical)
	freqDict = computeFreqDict(inDir, canInd)
	setWeightMatrices('PAM250', 'PAM250')
	topk = 25
	canonAnton = {1: [3], 2: [2,3], 3: [0,1]}
	triples = {1: [1,2,3], 2: [0,1,2], 3: [0,1,2]}
	singles = {1: [3], 2: [2], 3: [0]}
	
	canProt = 'RDLR'
	print canProt
	nmat = lookupCanonZF(freqDict, canProt, useNN = False, skipExact = False, 
	                     decompose = None, topk = None)
		
	print "Lookup:"
	print getConsensus(nmat)
	print "Final Matrix:"
	print nmat
	
	for k in [15, 20, 25, 30, 35, 40]:
		if k != 15:
			continue

		nmat, npb = lookupCanonZF(freqDict, canProt, useNN = True, skipExact = True, 
	                     	      decompose = singles, topk = k)
		
		print "Top %d: " %k
		print getConsensus(nmat)
		print "Final Matrix:"
		print nmat
				                 	

if __name__ == '__main__':
	main()
