import os
import re
import math
from jellyfish import hamming_distance

def getTripletPairEnrichment(tripDict, bgFreqs):
	# For each pair of triplets, compare the frequency
	# of variation in position -1, 2, 3, and 6 for all 
	# helices HD-1 from each other against
	# the background variation for these same positions.
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
	# lg( (freq of -1 variation for triplet-pair AAA,TTT)/
	#     (freq of -1 variation across the dataset) )

	tripPairDict = {}
	pMap = {0: -1, 1: 2, 2: 3, 3: 6}
	for i, t1 in enumerate(tripDict.keys()):
		#print i
		for t2 in tripDict.keys()[i+1:]:

			freqDict = {-1: [0, 0], 2: [0, 0], 
									 3: [0, 0], 6: [0, 0], 'numHD1Pairs': 0}
			totWeight = 0.0
			numHD1Pairs = 0
			for h1 in tripDict[t1]:
				for h2 in tripDict[t2]:
					if hamming_distance(h1, h2) == 1:
						for i, p1 in enumerate(h1):
							if p1 != h2[i]:
								freqDict['numHD1Pairs'] += 1
								freqDict[pMap[i]][1] += 1
								freqDict[pMap[i]][0] += tripDict[t1][h1] + tripDict[t2][h2]
								totWeight += tripDict[t1][h1] + tripDict[t2][h2]
								break
			
			# Normalize weights to frequencies and compute
			# the log ratios
			for k in sorted(freqDict.keys())[:-1]:
				try:
					freqDict[k][0] = freqDict[k][0]/totWeight
				except ZeroDivisionError:
					#print "No edges for %s %s (%d, %d)" \
					#	%(t1, t2, len(tripDict[t1]), len(tripDict[t2]))
					for k in sorted(freqDict.keys())[:-1]:
						print k, freqDict[k]
						freqDict[k][0] = 0
					break
				try:
					freqDict[k][0] = math.log(freqDict[k][0]/bgFreqs[k], 2)
				except ValueError:
					#print "No %d edges for %s %s (%d, %d)" \
					#	%(k, t1, t2, len(tripDict[t1]), len(tripDict[t2]))
					freqDict[k][0] = -100

			tripPairDict[(t1, t2)] = freqDict

	return tripPairDict

def getTripletHelices(dpath, trip, 
                      useWeights = False):
	# Map each helix to frequency across a triplet

	helixDict = {}
	fin = open(dpath + '4pos_' + trip + '.txt', 'r')
	for line in fin:
		#print(line)
		sp_line = line.strip().split()
		helix = sp_line[0]
		weight = eval(sp_line[1])
		if useWeights:
			helixDict[helix] = weight
		else:
			helixDict[helix] = 1
	fin.close()
	return helixDict

def getAllTripletHelices(dpath, triplets, 
                         useWeights = False):
	# Returns a dictionary of dictionaries, the outer
	# dictionary mapping each triplet to a helix
	# dictionary, which maps helices for that 
	# triplet to the associated frequencies.

	tripDict = {}
	for trip in triplets:
		#print trip
		tripDict[trip] = getTripletHelices(dpath, trip,
		                                   useWeights)
		#print(len(tripDict[trip]))
	#print(len(tripDict))
	return tripDict

def getHelixDict(dpath, useWeights = False):
	# map each helix to frequency across dataset

	helixDict = {}
	totWeight = 0.0
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
				if useWeights:
					helixDict[helix] += weight
					totWeight += weight
				else:
					pass
			else:
				if useWeights:
					helixDict[helix] = weight
					totWeight += weight
				else:
					helixDict[helix] = 1
		fin.close()

	# Normalize weights to a distribution
	for k in helixDict.keys():
		if useWeights:
			helixDict[k] = helixDict[k]/totWeight
		else:
			helixDict[k] = helixDict[k]/float(len(helixDict))

	return helixDict
	


def getHD1BGFreqs(dpath, useWeights = False):
	# Find the background percent of substituitions
	# that occur in position -1, 2, 3, and 6 among 
	# all canonical helices that are exactly hamming 
	# distance 1 apart across the entire dataset.  
	# If useWeights, then dataset-wide frquencies 
	# for each helix-pair when computing these bg frequencies.
	#
	# Return type is a dictionary mapping helix postion
	# to frequency

	helixDict = getHelixDict(dpath, useWeights)
	
	pMap = {0: -1, 1: 2, 2: 3, 3: 6}
	freqDict = {-1: 0, 2: 0, 3: 0, 6: 0}
	totWeight = 0.0
	for a, k1 in enumerate(helixDict.keys()):
		for b, k2 in enumerate(helixDict.keys()):
			if b < a:
				continue
			if hamming_distance(k1, k2) == 1:
				for i, p1 in enumerate(k1):
					if p1 != k2[i]:
						freqDict[pMap[i]] += helixDict[k1] + helixDict[k2]
						totWeight += helixDict[k1] + helixDict[k2]
						break

	# Normalize the bg frequencies
	for k in freqDict.keys():
		freqDict[k] = freqDict[k]/totWeight

	return freqDict

def makeTripDictTables(outPath, fing, strin, filt, 
                       tripPairDict, useWeights):
	
	if useWeights:
		wtag = 'weights'
	else:
		wtag = "noWeights"

	for pos in [-1, 2, 3, 6]:
		fout = open(outPath + '_'.join([fing, strin,
		            filt, wtag, 'a' + str(pos)]), 'w')
		header = "%s\t%s\t%s\t%s\t%s\n" %('t1', 't2', 'score',
		                                  'HD1pos', 'HD1allpos')
		fout.write(header)
		for (t1, t2) in sorted(tripPairDict.keys()):
			outStr = "%s\t%s\t%f\t%d\t%d\n" %(t1, t2,
			                          tripPairDict[(t1,t2)][pos][0],
			                          tripPairDict[(t1,t2)][pos][1],
			                          tripPairDict[(t1,t2)]['numHD1Pairs'])
			fout.write(outStr)
		fout.close()

def main():
	fing = 'F2'
	strin = 'union'
	filt = 'filt_10e-4_025_0_c'
	useWeights = False
	outPath = '../figures/positionVariation/'
	dpath = '../data/b1hData/antonProcessed/'
	dpath = '/'.join([dpath, fing, strin, filt])+'/'
	
	bases = ['A', 'C', 'G', 'T']
	triplets = []
	for b1 in bases:
		for b2 in bases:
			for b3 in bases:
				triplets.append(b1+b2+b3)

	bgFreqs = getHD1BGFreqs(dpath, useWeights = useWeights)
	print bgFreqs
	tripDict = getAllTripletHelices(dpath, triplets, 
	                                useWeights = useWeights)
	#print(len(tripDict.keys()))
	#for k in sorted(tripDict.keys()):
	#	print k, len(tripDict[k])
	
	tripPairDict = getTripletPairEnrichment(tripDict, bgFreqs)
	makeTripDictTables(outPath, fing, strin, filt, tripPairDict,
	                   useWeights = useWeights)

	#print bgFreqs
	#print sorted(tripPairDict.keys())
	#print len(tripPairDict)
	#for pair in tripPairDict.keys():
	#	summa = 0.0
	#	for k in tripPairDict[pair]:
	#		summa += tripPairDict[pair][k]
	#	print len(tripPairDict[pair]), summa


if __name__ == '__main__':
	main()