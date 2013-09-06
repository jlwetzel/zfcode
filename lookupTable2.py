import os
import re
import numpy as np
import math
from pwm import makeLogo, pwmfile2matrix, comparePCC, getConsensus, makeNucMatFile
from fixTables import normalizeFreq
from entropy import *
#from gatherBindStats import getProtDict

def getSubDict(fname, expo = False):
	# Return a substitution dictionary indicated
	# by the given file path.

	subDict = {}

	fin = open(fname, 'r')
	xlabs = fin.readline().strip().split()
	for line in fin:
		ylab = line.strip().split()[0]
		scores = [eval(i) for i in line.strip().split()[1:]]
		for i, s in enumerate(scores):
			if expo:
				subDict[ylab, xlabs[i]] = math.exp(s)
			else:
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
					nScore = weightMat[prot[i], a]
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

def sortAndWeightNeighbors(prot, bpos, neighbors, decomp, skipExact):
	# Sorts the neighbors in the order in which we 
	# would like to process them (most "important"
	# neighbors first).
	# 
	# Returns the newly sorted list of neighbors and 
	# a list of weights in the same order.
	
	# The order in which aminos should be allowed to vary 
	# for each of the bases.
	if len(decomp[1]) == 1:     # singles decomposition
		order = {1: [0, 1, 2],  # -1,2,3
			 	 2: [1, 0, 3],  #  2,-1,6
			 	 3: [3, 2, 1]}  #  6,3,2
	elif len(decomp[1]) == 2:   # doubles decomposition
		order = {1: [0, 1],     # -1,2
			 	 2: [1, 0],     #  2,-1
			 	 3: [3, 2]}     #  6,3
	elif len(decomp[1]) == 3:   # triples decomposition
		order = {1: [0],        # -1
			 	 2: [1],        #  2
			 	 3: [3]}        #  6

	
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
		orderInd = [(w,i) for (i,w) in enumerate(order_weights)]
		orderInd.sort(reverse = True)
		orderInd = [i for (w,i) in orderInd]
		
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

def getTopKNeighborsPWM(freqDict, prot, neighborDict, decomp, topk,
                        skipExact, verbose = None):
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
		
		# Sort the neighbors in the order in which we'd like 
		# to use them.
		if topk != None:
			#print k
			neighborDict[k], nWeights = sortAndWeightNeighbors(prot, k,
			                                                   neighborDict[k],
			                                                   decomp,
			                                                   skipExact)
		else:
			if NEIGHBOR_WEIGHTS != None:
				nWeights = getNeighborWeights(prot, neighborDict[k], skipExact)
			else:
				unifWeight = 1/float(len(neighborDict[k]))
				nWeights = np.array([unifWeight]*len(neighborDict[k]),
			                       dtype = float)

		#print k
		#print neighborDict[k]
		# Get the normalized frequency vectors for each neighbor at  
		# this base position
		#print
		#print len(neighborDict[k])
		for i, n in enumerate(neighborDict[k]):
			newVect = getNeighborBaseVect(freqDict, n, k)
			baseVectors.append(newVect)
		
		# Remove Nones from both the neighbor list, baseVectorList, and the 
		# corresponding entries in nWeights, then renormalize nWeights
		nWeights = list(nWeights)
		nWeights = np.array([nWeights[i] for i, n \
		                    in enumerate(baseVectors) if n != None][:topk], 
		                    dtype = float)
		nWeights = nWeights/np.sum(nWeights)
		neighborDict[k] = [neighborDict[k][i] for i, n \
		                    in enumerate(baseVectors) if n != None][:topk]
		baseVectors = [i for i in baseVectors if i != None][:topk]
		neighborsPerBase.append(len(baseVectors))

		# If we found at least one neighbor
		if baseVectors != []:
			# Add the weighted vectors together
			if verbose != None:
				print
			for i in range(len(baseVectors)):
				pwm[k-1] = pwm[k-1] + baseVectors[i]*nWeights[i]		
				
				if verbose != None:
					print "Base: %d\tNeighbor: %s\tWeight: %.5f" \
					%(k, neighborDict[k][i], nWeights[i])
					print baseVectors[i]
					dirname = verbose
					searchStr = '"%s.%s%s.%s"'%(neighborDict[k][i][0],
					                            neighborDict[k][i][1],
					                            neighborDict[k][i][2],
					                            neighborDict[k][i][3])
					handle = os.popen("grep %s %s*.txt" %(searchStr, dirname))
					for line in handle:
						sp_line = line.strip().split(':')
						file = sp_line[0].split('/')[-1]
						line_info = '  '.join(sp_line[1].split('\t')[:4])
						if file == 'all.txt':
							continue
						print "%s: %s" %(file, line_info)

			# Renormalize since some neighbors may not have been used
			pwm[k-1] = pwm[k-1]/np.sum(pwm[k-1]) 
			if verbose != None:
				print "Used %d neighbors for base %d" \
					%(neighborsPerBase[len(neighborsPerBase)-1], k)

		
		# If we found no neighbors use a uniform vector
		else:
			pwm[k-1] = np.array([0.25, 0.25, 0.25, 0.25], dtype = float)

	#print neighborsPerBase
	return pwm, neighborsPerBase


def get3merList(freqDict, protein, canonical = False,
				useNN = False, skipExact = False, 
				decompose = None, topk = None, 
				verbose = None):
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

	return getTopKNeighborsPWM(freqDict, protein, neighborDict, decompose,
	                           topk, skipExact, verbose)
	
def makeDir(path):
	# Try to make a path and pass if it can't be created
	try:
		os.mkdir(path)
	except OSError:
		pass

def lookupCanonZFArray(freqDict, canonZFs, useNN = True, 
                       skipExact = False, decompose = None, topk = None,
                       verbose = None):
	# Performs modular lookup for an array of canonical 
	# helix-position ZF domains.  Domains should be given
	# in reverse order of their appearance on the protein.

	numZFs = len(canonZFs)
	pwm = np.zeros(shape = (numZFs*3, 4))
	
	for i in range(numZFs):
		nmat = lookupCanonZF(freqDict, canonZFs[i], 
		                     useNN, skipExact, decompose, topk, verbose)
		for j in range(len(nmat)):
			pwm[i*len(nmat) + j,:] = nmat[j,:]
	
	return pwm

def lookupCanonZF(freqDict, canonProt, useNN = True, skipExact = False,
                  decompose = None, topk = None, verbose = None):
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
	                                         useNN, skipExact, decompose, topk, 
	                                         verbose)
	
	# This is target/frequncy list
	if isinstance(targList, list):
		# If the target list is empty, we can't return a specificity,
		# so we just return a uniform distribution for each base
		if targList == []:
			nucmat = np.zeros((3,4), float)
			for i in range(len(nucmat)):
				for j in range(len(nucmat[i,:])):
					nucmat[i,j] = 0.25
			return nucmat, [0,0,0]

		# Convert the target list to a position freq mat and return
		return targListToFreqMat(targList), [0,0,0]

	# This is already a proper numpy array
	else:
		return targList, neighborsPerBase

def setWeightMatrices(weight_mat, order_mat):
	# Set the global weight matrix parameters
	global NEIGHBOR_WEIGHTS
	global NEIGHBOR_ORDER
	
	matLoc = '../data/substitution_mats/'

	print weight_mat.split('.')[0].split('_')[-1]
	if weight_mat.split('.')[0].split('_')[-1] == 'bin':
		expoW = False
	else:
		expoW = True

	if order_mat.split('.')[0].split('_')[-1] == 'bin':
		expoO = False
	else:
		expoO = True

	NEIGHBOR_WEIGHTS = getSubDict(matLoc + weight_mat + '.txt', 
	                              expo = expoW)
	NEIGHBOR_ORDER = getSubDict(matLoc + order_mat + '.txt',
	                            expo = expoO)	

def getPosEntropies(freqDict, norm = False):
	# Working on this

	entropyDict = {}
	napos = len(freqDict[freqDict.keys()[0]].keys()[0])
	aseqFreqDict = {}
	for nseq in freqDict.keys():
		for aseq in freqDict[nseq].keys():
			if aseqFreqDict.has_key(aseq):
				aseqFreqDict[aseq] += freqDict[nseq][aseq]
			else:
				aseqFreqDict[aseq] = freqDict[nseq][aseq]
	for i in range(napos):
		ifreqDict = {}
		for aseq in aseqFreqDict.keys():
			if ifreqDict.has_key(aseq[i]):
				ifreqDict[aseq[i]] += aseqFreqDict[aseq]
			else:
				ifreqDict[aseq[i]] = aseqFreqDict[aseq]
		freqs = np.array([f for (k, f) in aseqFreqDict.items()], dtype = float)
		freqs = freqs/np.sum(freqs)
		entropyDict['a' + repr(i)] = entropy(freqs)

	print entropyDict
	

def main():

	# Debugging stuff
	outDir = '../data/scratch/'
	inDir = '../data/b1hData/newDatabase/6varpos/F2/low/protein_seq_cut10bc_0_5/'
	canonical = True
	varpos = 6
	canInd = getPosIndex(varpos, canonical)
	freqDict = computeFreqDict(inDir, canInd)
	#setWeightMatrices('PAM250', 'PAM250')
	setWeightMatrices('PAM30', 'PAM30')
	topk = 25
	
	triples = {1: [1,2,3], 2: [0,2,3], 3: [0,1,2]}
	doubles = {1: [2,3], 2: [2,3], 3: [0,1]}
	singles = {1: [3], 2: [2], 3: [0]}
	decomp = singles
	
	# Working on this
	#entropyDict = getPosEntropies(freqDict)
	f3s = [('ATT','FQSGLIQ'), ('CAC', 'HQANLIH'),('CGG', 'RNAHLTD'),
		   ('GAC', 'DQSNLTR'), ('GCG', 'TKYDLTR'), ('TGT', 'MKQHLTY'),
		   ('AAA', 'SAGSLYN'), ('CAA', 'QKINLIN'), ('ATC', 'DKSYLYT'),
		   ('CCC', 'LKKTLTD'), ('CTC', 'DQSALLG'), ('GCC', 'LKKTLIN')]
	
	f2s = [('ACC', 'CPKALRA'), ('CCC', 'DPRCLVR'), ('CGA', 'SYNGLRG'),
		   ('GTT', 'FAMQLTR'), ('TCC', 'DPRSLVY'), ('TAA', 'SYMGLRG')]
	hughes = ['DCRDLAR', 'TSGELVR', 'TSGHLVR', 'QSSHLTR', 'RSDDLQR', 'QSGHLQR',
			  'TSGNLVR', 'QKSSLIA', 'TTGNLTV', 'SPADLTR', 'QKSSLIA', 'TSGELVR']
	
	hughesTargs = ['XXX'] * len(hughes)
	hughes = zip(hughesTargs, hughes)
	scratch = [('XXX', 'RXDLXXR')]


	pairs = scratch
	targs = [i[0] for i in pairs]
	prots = [i[1] for i in pairs]
	canProts = [i[0]+i[2]+i[3]+i[6] for i in prots]


	for i, canProt in enumerate(canProts):
		print targs[i], canProt
		nmat, npb = lookupCanonZF(freqDict, canProt, useNN = False, skipExact = False, 
		                     decompose = None, topk = None, verbose = None)
		
		label = '_'.join([targs[i], prots[i], 'lookup'])
		makeNucMatFile(outDir, label, nmat)
		logoIn = outDir + label + '.txt'
		logoOut = outDir + label + '.pdf'
		makeLogo(logoIn, logoOut, alpha = 'dna', 
			         colScheme = 'classic',
			         annot = "'5,M,3'",
			         xlab = '_'.join([targs[i], prots[i]]))
		
		print "Lookup:"
		print getConsensus(nmat)
		print "Final Matrix:"
		print nmat
		
		for k in [10, 15, 20, 25, 30, 35, 40]:


			if k != 10:
				continue

			nmat, npb = lookupCanonZF(freqDict, canProt, useNN = True, skipExact = False, 
		                     	      decompose = decomp, topk = k, verbose = inDir)
			
			label = '_'.join([targs[i], prots[i], 'top' + str(k)])
			makeNucMatFile(outDir, label, nmat)
			logoIn = outDir + label + '.txt'
			logoOut = outDir + label + '.pdf'
			makeLogo(logoIn, logoOut, alpha = 'dna', 
				         colScheme = 'classic',
				         annot = "'5,M,3'",
				         xlab = '_'.join([targs[i], prots[i]]))

			print "Top %d: " %k
			print getConsensus(nmat)
			print "Final Matrix:"
			print nmat	                 	

if __name__ == '__main__':
	main()
