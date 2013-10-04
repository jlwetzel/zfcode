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
		freq = eval(sp_line[1])
		rest = '\t'.join(sp_line[2:])
		seqDict[prot] = freq

	return seqDict

def unionSeqDicts(d1, d2):
	# Takes as input a pair of seqdicts and unions them
	# averaging the frequencies and outputting the 
	# rest of the line for the first sequence
	#
	# Each dict is of for seq -> freq
	# Returns a dictionary of the same form

	unionDict = {}
	unionKeys = set(d1.keys()) | set(d2.keys())
	for k in unionKeys:
		if d1.has_key(k) and d2.has_key(k):
			unionDict[k] = (d1[k] + d2[k])/2
		elif d1.has_key(k):
			unionDict[k] = d1[k]/2
		else:
			unionDict[k] = d2[k]/2
	return unionDict

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


def combineHighAndLow(path, filt):
	# Combines the high and low stringency datasets,
	# taking a union of sequences and averaging.
	# Path is to the finger directory
	# filter is which filter to be combined on 
	# high and low stringency

	hpath = path + 'high/' + filt + '/'
	lpath = path + 'low/' + filt + '/'
	# Make new directories for combined files
	newPath = path + '/union/' + filt + '/'
	makeDir(path + '/union/')
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
		unionDict = unionSeqDicts(highDict, lowDict)

		summa = 0
		for k in unionDict.keys():
			summa += unionDict[k]
		print summa

		# Write to a file
		writeToFile(unionDict, newPath + fname)

def main():

	"""
	path = '../data/b1hData/antonProcessed/F2/'
	filt = 'filt_10e-4_025_0_c'
	combineHighAndLow(path, filt)
	"""

	fings = ["F1", "F2", "F3"]
	filts = ["filt_10e-4_025_0_c"]
	for f in fings:
		for filt in filts:
			path = '../data/b1hData/antonProcessed/' + f + '/'
			combineHighAndLow(path, filt)

if __name__ == '__main__':
	main()