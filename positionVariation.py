import os
import re
import math
from jellyfish import hamming_distance

def getTripletPairEnrichment(tripDict, bgFreqs,
                             weightObs = False):
	# For each pair of triplets, compare the weight 
	# of observed HD-1 variations in postions 1, 2,
	# 3 or 6 to the weight of expected variations.
	#
	# tripDict is a dictionary of dictionaries, 
	# mapping triplets to dicts which in turn map 
	# helices to their per-triplet frequencies.
	# bgFreq is also a dict, mapping -1, 2, 3, and 6 to
	# the corresponding bg freqs.
	#
	# Return value is a dictionary of dictionaries,
	# the outer dictionary mapping tuple-pairs of triplets
	# to dictionaries which map positions (-1,2,3, or 6) to
	# enrichment scores:  e.g.:
	# lg( (obs. weight of -1 variation for triplet-pair AAA,TTT)/
	#     (exp. weight of -1 variation for triplet-pair AAA,TTT) )

	tripPairDict = {}
	obsDict = {}
	expDict = {}
	pMap = {0: -1, 1: 2, 2: 3, 3: 6}
	for a, t1 in enumerate(sorted(tripDict.keys())):
		for b, t2 in enumerate(sorted(tripDict.keys())):

			# Don't calculate above diagonal
			if b < a:
				tripPairDict[(t1,t2)] = {-1: 0, 2: 0, 3: 0, 6: 0}
				obsDict[(t1,t2)] = {-1: -1, 2: -1, 3: -1, 6: -1}
				expDict[(t1,t2)] = {-1: -1, 2: -1, 3: -1, 6: -1}
				continue

			# Get weight of observed variations, total 
			# weight for all pairs, and minimum weight 
			# of any possible pair
			obsPosWeights = {-1: 0, 2: 0, 3: 0, 6: 0}
			totWeight = 0.0
			minWeight = 1.5
			for h1 in tripDict[t1]:
				for h2 in tripDict[t2]:
					hPairWeight = tripDict[t1][h1] * tripDict[t2][h2]
					if hamming_distance(h1, h2) == 1:
						for i, p1 in enumerate(h1):
							if p1 != h2[i]:
								if weightObs:
									obsPosWeights[pMap[i]] += hPairWeight
								else:
									obsPosWeights[pMap[i]] += 1
								break
					if weightObs and (hPairWeight < minWeight):
						minWeight = hPairWeight
					totWeight += hPairWeight
			
			# Do "plus 1" smoothing
			num0Freqs = 0
			for k in obsPosWeights.keys():
				if obsPosWeights[k] == 0:
					num0Freqs += 1
			if num0Freqs > 0:
				for k in obsPosWeights.keys():
					if weightObs:
						obsPosWeights[k] += minWeight
					else:
						obsPosWeights[k] += 1

			# Get number of the expected variations of each type
			expDict[(t1,t2)] = {}
			for k in obsPosWeights.keys():
				expDict[(t1,t2)][k] = bgFreqs[k]*totWeight

			# Compare observed to expected
			tripPairDict[(t1,t2)] = {}
			for k in obsPosWeights.keys():
				tripPairDict[(t1,t2)][k] = \
					math.log((obsPosWeights[k]/ \
					         expDict[(t1,t2)][k]), 2)

			obsDict[(t1,t2)] = obsPosWeights

	return tripPairDict, obsDict, expDict

def getTripletHelices(dpath, trip, 
                      helixWeights):
	# Map each helix to frequency across a triplet

	helixDict = {}
	fin = open(dpath + '4pos_' + trip + '.txt', 'r')
	for line in fin:
		#print(line)
		sp_line = line.strip().split()
		helix = sp_line[0]
		weight = eval(sp_line[1])
		helixDict[helix] = weight/(helixWeights[helix])
	fin.close()
	return helixDict

def getAllTripletHelices(dpath, triplets, helixWeights):
	# Returns a dictionary of dictionaries, the outer
	# dictionary mapping each triplet to a helix
	# dictionary, which maps helices for that 
	# triplet to the associated frequencies.

	tripDict = {}
	for trip in triplets:
		tripDict[trip] = getTripletHelices(dpath, trip,
		                                   helixWeights)
	return tripDict	

def getHelixDict(dpath):
	# Get the set of helices in across the entie datset
	# along with their dataset-wide weights

	helixDict = {}
	handle = os.popen('ls %s' %dpath)
	for fname in handle:
		fname = fname.strip()
		if re.match(r'^4pos_[ACGT]{3}.txt', fname) == None:
			continue

		fin = open(dpath + fname, 'r')
		for line in fin:
			sp_line = line.strip().split()
			helix = sp_line[0]
			weight = eval(sp_line[1])
			if helixDict.has_key(helix):
				helixDict[helix] += weight
			else:
				helixDict[helix] = weight
		fin.close()

	return helixDict

def getHD1BGProbs(dpath):
	# Estimate the probability that a pair of helices
	# have hamming distance 1 and differ in the -1, 2, 
	# 3, or 6 positions by taking finding the frequency
	# with which this happens in our dataset

	helixDict = getHelixDict(dpath)
	helixList = helixDict.keys()

	pMap = {0: -1, 1: 2, 2: 3, 3: 6}
	hd1Counts = {-1: 0, 2: 0, 3: 0, 6: 0}
	hd1Freqs = {-1: 0, 2: 0, 3: 0, 6: 0}
	for a, k1 in enumerate(helixList):
		for k2 in helixList[a+1:]:
			if hamming_distance(k1, k2) == 1:
				for i, p1 in enumerate(k1):
					if p1 != k2[i]:
						hd1Counts[pMap[i]] += 1
						break

	# Normalize the bg frequencies
	numPosEdges = (len(helixList) * (len(helixList) - 1))/2
	for k in hd1Counts.keys():
		hd1Freqs[k] = float(hd1Counts[k])/numPosEdges

	return hd1Freqs, helixDict

def makeTripDictTables(outPath, fing, strin, filt, 
                       tripPairDict, obsDict, expDict,
                       weightObs = False):

	if weightObs:
		wtag = "weightObs"
	else:
		wtag = "noWeightObs"

	for pos in [-1, 2, 3, 6]:
		fout = open(outPath + '_'.join([fing, strin,
		            filt, wtag, 'a' + str(pos)]) + '.txt', 'w')
		header = "%s\t%s\t%s\t%s\t%s\n" %('t1', 't2', 'score',
		                                  'obs', 'exp')
		fout.write(header)
		for (t1, t2) in sorted(tripPairDict.keys()):
			outStr = "%s\t%s\t%f\t%f\t%f\n" \
				%(t1, t2, tripPairDict[(t1,t2)][pos], \
				  obsDict[(t1,t2)][pos],expDict[(t1,t2)][pos])
			fout.write(outStr)
		fout.close()

def main():
	fing = 'F2'
	strin = 'union'
	filt = 'filt_10e-4_025_0_c'
	outPath = '../figures/positionVariation/'
	dpath = '../data/b1hData/antonProcessed/'
	dpath = '/'.join([dpath, fing, strin, filt])+'/'
	weightObs = True
	
	bases = ['A', 'C', 'G', 'T']
	triplets = []
	for b1 in bases:
		for b2 in bases:
			for b3 in bases:
				triplets.append(b1+b2+b3)

	bgProbs, helixWeights = getHD1BGProbs(dpath)
	print bgProbs
	tripDict = getAllTripletHelices(dpath, triplets,
	                                helixWeights)
	tripPairDict, obsDict, expDict = \
		getTripletPairEnrichment(tripDict, bgProbs, 
		                         weightObs = weightObs)
	makeTripDictTables(outPath, fing, strin, filt, 
	                   tripPairDict, obsDict, expDict,
	                   weightObs = weightObs)


if __name__ == '__main__':
	main()