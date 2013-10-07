# Creates vectors of binding probabilities across the 64 
# possible target 3-mers for each of the 
# filtered datasets unique canonical binders and writes them
# to files in the filter's directory.
# Also computes the entropy for each of these distributions
# and writes them to the same file

from itertools import product
from gatherBindStats import getProtSet
from lookupTable2 import computeFreqDict, get3merList
import numpy as np

def computeEntropy(f):
    # Computes a normalized Shannon Entropy for the 
    # distribution f in base 2.  Follows the convention that 
    # plog(p) = 0 if p = 0.

    if len(f) == 1:
        return 0.0

    fentropy = 0.0
    for i in range(len(f)):
        if f[i] == 0:
            continue
        else:
            fentropy -= f[i] * np.log2(f[i])

    return fentropy

def getProtWeights(path, canInd):
	# Maps each unique protein to its weight across the 
	# entire dataset as well as its number of coding variants	
	# path should be to a file of the format 'all.txt'
	# canInd should be to the set of indicies which 
	# we are interrested in

	fin = open(path, 'r')
	protWeightDict = {}
	for line in fin:
		sp_line = line.strip().split()
		targ, prot = sp_line[0], sp_line[1]
		freq = eval(sp_line[2])
		canProt = ''
		for i in canInd:
			canProt += prot[i]

		if protWeightDict.has_key(canProt):
			protWeightDict[canProt] += freq
		else:
			protWeightDict[canProt] = freq

	# Normalize to sum to 1
	for k in protWeightDict.keys():
		protWeightDict[k] /= 64.0

	return protWeightDict

def getBindVects(protSet, bindDir, canInd, allTargs):
	# Returns a dictionary of normalized binding frequency
	# list for each prot in protSet.  The prots are the 
	# keys and the vectors are the values.  
	# allTargs the ordered list of all possible targets
	
	freqDict = computeFreqDict(bindDir, canInd)
	vectDict = {}
	for prot in protSet:
		targList, npb = get3merList(freqDict, prot, canonical = True)
		targsHit = [i[0] for i in targList]
		targsFreq = [i[1] for i in targList]
		#print targsHit
		#print targsFreq
		vectDict[prot] = []
		for t in allTargs:
			try:
				tInd = targsHit.index(t)
				vectDict[prot].append(targsFreq[tInd])
			except ValueError:
				vectDict[prot].append(0.0)
		vectDict[prot] = np.array(vectDict[prot], dtype = 'float')

	return vectDict

def makeBindVectFiles(bindPref, fings, strins, filts, canInd, allTargs):
	# Compute the normalized vectors for binding across all 
	# of the 64 possible targets and write to file in the binding 
	# directory

	for f in fings:
		for s in strins:
			for filt in filts:
				
				# Get all of the binding vectors
				bindDir = bindPref + '/' + f + '/' + s + '/' + \
					filt + '/'
				fpath = bindDir + 'all.txt'
				protSet = getProtSet(fpath, canInd)
				protWeights = getProtWeights(fpath, canInd)
				bindVects = getBindVects(protSet, bindDir, canInd, allTargs)
				
				# Open a file for this binding directory and write header
				fout = open(bindDir + 'bindVectsUnique4pos.txt', 'w')
				head = "prot\tentropy\tweight"
				for targ in allTargs:
					head += '\t' + targ
				head += '\n'
				fout.write(head)

				# Get entropy and complete line to file
				weightedKeys = [(protWeights[k], k) for k in bindVects.keys()]
				weightedKeys.sort(reverse = True)

				for prot in [i[1] for i in weightedKeys]:					
					print prot
					# Print the line to file
					outStr = prot
					entropy = computeEntropy(bindVects[prot])
					outStr += '\t' + str(entropy) + '\t' + str(protWeights[prot])
					for i in range(len(bindVects[prot])):
						outStr += '\t' + str(bindVects[prot][i])
					outStr += '\n'
					fout.write(outStr)

				# Close the file
				fout.close()

def main():
	bindPref = '../data/b1hData/antonProcessed/'
	fings = ["F1", "F2", "F3"]
	strins = ["inter"]
	filts = ['filt_10e-4_025_0_c']
	canInd = [0,2,3,5]

	# Make a list of all triplets
	bases = ['A', 'C', 'G', 'T']
	allTargs = []
	for i in product(bases, bases, bases):
		allTargs.append(''.join(i))

	# Make binding vector files
	makeBindVectFiles(bindPref, fings, strins, filts, canInd, allTargs)

if __name__ == '__main__':
	main()

