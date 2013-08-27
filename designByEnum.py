from lookupTable2 import computeFreqDict, getPosIndex, lookupCanonZF
import os, re
from pwm import makeNucMatFile, makeLogo

def makeDir(path):
	# Try to make a path and pass if it can't be created
	try:
		os.mkdir(path)
	except OSError:
		pass

def getTopProts(inDir, topkPer3mer, ind):
	# Gathers the top k best binding unique proteins with 
	# per 3mer with uniqueness determinied by the indices ind.
	# topkPer3mer is a list that gives the value of k for 
	# each of the 3-mers.
	# Returns the proteins as a dictionary of ordered lists, 
	# using the 3-mer as the dictionary key.

	handle = os.popen('ls ' + inDir)

	# For each 3-mer file in this directory

	topProtDict = {}
	targNum = 0
	for fname in handle:
		fname = fname.strip()

		# Skip files that are not named properly
		if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
			continue

		targ = fname.split('.')[0]  # The 3-mer target for this file
		orderedList = []  # List of proteins to be return for this 3mer
						  # in decreasing order of max frequency
		protSeqs = {}  # Just a hash table for quickly checking if we have a
					   # protein sequence already  
		
		fin = open(inDir + fname, 'r')
		numProtsFound = 0   # Number of unique prots found so far for this 3mer
		for lNum, line in enumerate(fin):

			if numProtsFound == topkPer3mer[targNum]:
				targNum += 1
				break

			# Get protein with all varied positions and frequency
			sp_line = line.strip().split()
			wholeProt = sp_line[0]
			freq = eval(sp_line[1])

			# Get the protein for the indices supplied
			prot = ''
			for i in ind:
				prot += wholeProt[i]

			# If this is a new protein then record and append to orderd list
			if protSeqs.has_key(prot):
				pass
			else:
				orderedList.append([prot, freq, lNum])
				protSeqs[prot] = numProtsFound
				numProtsFound += 1

		fin.close()
		topProtDict[targ] = orderedList

	return topProtDict

def predictTopProts(freqDict, topProtDict, outDir, topk = 25):
	# Makes DNA binding predictions for each of the proteins
	# in topProtDict using the given frequency dictionary
	# and places the pwms and logos in outDir in an organized way.

	singles = {1: [3], 2: [2], 3: [0]}

	for targ in sorted(topProtDict.keys()):
		orderedProts = [i[0] for i in topProtDict[targ]]
		orderedFreqs = [i[1] for i in topProtDict[targ]]
		orderedLines = [i[2] for i in topProtDict[targ]]
		
		targDir = outDir + targ + '/'
		print targDir
		makeDir(targDir)
		pwmdir = targDir + 'pwms/'
		logodir = targDir + 'logos/'
		makeDir(pwmdir)
		makeDir(logodir)


		for i, prot in enumerate(orderedProts):
			
			# Make the nearest neighbors lookup prediction
			nmat = lookupCanonZF(freqDict, prot, useNN = True,
			                     skipExact = False, decompose = singles,
			                     topk = topk)
			
			# Create the frequency matrix file
			label = '_'.join([targ, prot, str(orderedLines[i])])
			makeNucMatFile(pwmdir, label, nmat)

			# Create the logo file
			logoIn = pwmdir + label + '.txt'
			logoOut = logodir + label + '.pdf'
			makeLogo(logoIn, logoOut, alpha = 'dna', 
		         colScheme = 'classic',
		         annot = "'5,M,3'",
		         xlab = '_'.join([targ,prot]))


def main():
	
	canonical = True
	varpos = 6
	canInd = getPosIndex(varpos, canonical)
	inDir = '../data/b1hData/newDatabase/6varpos/F2/low/protein_seq_cut10bc_0_5/'
	outDir = '../data/design/F2_low_nn25_cut10bc_0_5_topFrac_50/'
	topFrac = .50
	maxPerTarg = 100

	# Get the frequency dictionary for doing predictions
	freqDict = computeFreqDict(inDir, canInd)

	# Figure out how many proteins for top x% for each 3-mer
	topkPer3mer = [int(topFrac*(len(freqDict[k]))) \
		for k in sorted(freqDict.keys())]
	# Place a cap on how many we can look at
	topkPer3mer = [min(maxPerTarg, i) for i in topkPer3mer]

	# Get the top x percent of proteins with regard to
	# unique canoncial sequence
	topProtDict = getTopProts(inDir, topkPer3mer, canInd)
	for k in sorted(topProtDict.keys()):
		print k, len(topProtDict[k])

	# Create logos/pfms for each of the prots in topProtDict
	#predictTopProts(freqDict, topProtDict, outDir)

	#print freqDict.keys()
	#topkProtDict = geTopkProtDict(inDir)

if __name__ == '__main__':
	main()