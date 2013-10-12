import os
import re

def makeDir(path):
	# Try to make a path and pass if it can't be created
	try:
		os.mkdir(path)
	except OSError:
		pass

def getSeqDict(fname):
	# Returns a dictionary of sequence info lines 
	# for the data in the handle

	seqDict = {}
	fin = open(fname, 'r')
	for line in fin:
		sp_line = line.strip().split()
		prot = sp_line[0]
		freq = float(sp_line[1])
		rest = '\t'.join(sp_line[2:])
		seqDict[prot] = freq

	return seqDict

def combineSeqDicts(d1, d2, combineType):
	# Takes as input a pair of seqdicts and unions them
	# averaging the frequencies and outputting the 
	# rest of the line for the first sequence
	#
	# Each dict is of for seq -> freq
	# Returns a dictionary of the same form

	combineDict = {}

	if combineType == 'union':
		combineKeys = set(d1.keys()) | set(d2.keys())
	elif combineType == 'inter':
		combineKeys = set(d1.keys()) & set(d2.keys())

	for k in combineKeys:
		if d1.has_key(k) and d2.has_key(k):
			combineDict[k] = (d1[k] + d2[k])/2
		elif d1.has_key(k):
			combineDict[k] = d1[k]/2
		else:
			combineDict[k] = d2[k]/2
	
	return combineDict


def writeToFile(seqDict, path):
	# Writes a file of the form AAA.txt
	# seqDict is of the form targ -> freq
	# path is the location for the new file
	# Dummy values are inserted as placeholders
	# for numObs, numPoss, and entropy

	fout = open(path, 'w')

	# Convert to list and sort by frequency
	l = [[seqDict[k], k] for k in seqDict.keys()]
	l.sort(reverse = True)
	l = [[i[1], i[0]] for i in l]

	for i in l:
		fout.write("%s\t%f\t%d\t%d\t%d\n" \
		           	%(i[0], i[1], -1, -1, -1))

	fout.close()

def normalizeSeqDict(seqDict):
	# Takes a dictionary of seqs mapping to real 
	# values and normalizes the values so that 
	# they add to one.

	# Get the total
	tot = 0.0
	for k in seqDict.keys():
		tot += seqDict[k]

	# Normalize
	for k in seqDict.keys():
		seqDict[k] = seqDict[k]/tot

def combineHighAndLow(path, filt, combineType):
	# Combines the high and low stringency datasets,
	# taking a union of sequences and averaging.
	# Path is to the finger directory
	# filter is which filter to be combined on 
	# high and low stringency

	hpath = path + 'high/' + filt + '/'
	lpath = path + 'low/' + filt + '/'
	
	# Make new directories for combined files
	newPath = path + '/' + combineType + '/' + filt + '/'
	makeDir(path + '/' + combineType + '/')
	makeDir(newPath)

	highHandle = os.popen('ls %s' %(hpath))
	# Combine each pair of files (low and high) per target
	for fname in highHandle:
		fname = fname.strip()
		if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
			continue

		# Get the dictionaries for combining
		highDict = getSeqDict(hpath + fname)
		lowDict = getSeqDict(lpath + fname)

		# Combine by averaging frequencies
		combineDict = combineSeqDicts(highDict, lowDict,
		                              combineType)

		# Renormalize if we are dropping anything
		if combineType == 'inter':
			normalizeSeqDict(combineDict)

		# Write to a file
		writeToFile(combineDict, newPath + fname)

def combineFingers(path1, path2, filt, outPath, combineType,
                   interType = None):

	# Combines datasets across fingers
	# combineType is either 'union' or 'inter'
	# interType is only needed if combineType is 'inter'
	# and specifices whether to intersect at the level of 
	# triplet files or at the level of whole datasets

	path1 = path1 + '/' + filt + '/'
	path2 = path2 + '/' + filt + '/'
	outPath = outPath + '/' + filt + '/'
	makeDir(outPath)


	p1Handle = os.popen('ls %s' %(path1))
	
	# Handle the case of unioning or intersecting based on 
	# individual triplet files
	if (combineType == 'inter' and interType == "triples") \
		or combineType == 'union':

		# Combine each pair of files per target
		for fname in p1Handle:
			fname = fname.strip()
			if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
				continue

			# Get the dictionaries for combining
			p1Dict = getSeqDict(path1 + fname)
			p2Dict = getSeqDict(path2 + fname)

			# Combine by averaging frequencies
			combineDict = combineSeqDicts(p1Dict, p2Dict,
			                              combineType)

			# Renormalize if we are dropping anything
			if combineType == 'inter':
				normalizeSeqDict(combineDict)

			# Write to a file
			writeToFile(combineDict, outPath + fname)


	# Handle the case of intersecting on all unique proteins
	# in the two datasets.
	elif (combineType == 'inter' and interType == 'all'):

		# Get the intersection of unique canonical proteins
		# across the two fingers
		p1Set = set()
		p2Set = set()
		for fname in p1Handle:
			fname = fname.strip()
			if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
				continue

			p1Dict = getSeqDict(path1 + fname)
			p2Dict = getSeqDict(path2 + fname)
			for k in p1Dict.keys():
				p1Set.add(k)
			for k in p2Dict.keys():
				p2Set.add(k)

		print len(p1Set), len(p2Set)
		interSet = p1Set & p2Set
		print(len(interSet))

		# Combine each pair of files per target
		p1Handle = os.popen('ls %s' %(path1))
		for fname in p1Handle:
			fname = fname.strip()
			if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
				continue

			# Get the dictionaries for combining
			p1Dict = getSeqDict(path1 + fname)
			p2Dict = getSeqDict(path2 + fname)

			# Combine by averaging frequencies
			for k in p1Dict.keys():
				if k not in interSet:
					del p1Dict[k]
			for k in p2Dict.keys():
				if k not in interSet:
					del p2Dict[k]
			combineDict = combineSeqDicts(p1Dict, p2Dict,
			                              'union')

			# Renormalize if we are dropping anything
			if combineType == 'inter':
				normalizeSeqDict(combineDict)

			# Write to a file
			writeToFile(combineDict, outPath + fname)



def main():

	"""
	path = '../data/b1hData/antonProcessed/F2/'
	filt = 'filt_10e-4_025_0_c'
	combineHighAndLow(path, filt)
	"""

	"""
	# Combining high and low for each finger
	fings = ["F1", "F2", "F3"]
	filts = ["filt_10e-4_025_0_c"]
	for f in fings:
		for filt in filts:
			path = '../data/b1hData/antonProcessed/' + f + '/'
			combineHighAndLow(path, filt, 'inter')
	"""

	# Unioning between fingers F2 and F3 (intersections)
	filt = "filt_10e-4_025_0_c"
	path1 = '../data/b1hData/antonProcessed/F2/high/'
	path2 = '../data/b1hData/antonProcessed/F3/high/'
	outPath = '../data/b1hData/antonProcessed/F2F3/' + \
		'unionHigh/'
	combineFingers(path1, path2, filt, outPath,
	               combineType = 'union', interType = 'all')


if __name__ == '__main__':
	main()