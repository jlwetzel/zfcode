from lookupTable2 import computeFreqDict, getPosIndex, lookupCanonZF
import os, re
from pwm import makeNucMatFile, makeLogo, pwmfile2matrix
import numpy as np

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

def getTopTenProts(targ, nucMatDict, opt = 'avgMinDiff'):
	# nucMatDict is a dictionary of binding preference
	# frequency matrices, indexed by the filenames for
	# those frequency matrices.
	# targ is a 3mer target sequence
	#
	# This function returns a list of the top 10 filenames 
	# for optimizing the objective given by opt
	# (in decreasing order)

	nucs = ['A', 'C', 'G', 'T']
	top10 = {}
	fnames = [k for k in sorted(nucMatDict.keys())]
	#print fnames

	# Maximize the average (or minimim) minimum difference 
	# between the desired nucleotide's frequency and the 
	# next most frequent nucleotide's frequency 
	# across all the nucleotide positions. 
	if opt == 'avgMinDiff' or opt == 'minMinDiff':
		fnameScores = []
		for fname in fnames:
			nmat = nucMatDict[fname]
			
			minDiffs = []
			for i in range(len(nmat)):
				
				# Score diff as -1 if desired nuc doesn't have highest freq
				if np.argmax(nmat[i,:]) != nucs.index(targ[i]):
					minDiffs = [-1] * len(nmat)
					break

				targFreq = np.max(nmat[i,:])
				maxVal = 0.0
				for j in range(len(nmat[i,:])):
					if j != nucs.index(targ[i]) and nmat[i,j] > maxVal:
						maxVal = nmat[i,j]
				minDiffs.append(targFreq - maxVal)

			if opt == 'avgMinDiff':
				fnameScores.append(np.mean(minDiffs))
			elif opt == 'minMinDiff':
				fnameScores.append(min(minDiffs))


	#print fnameScores
	# Find the index order for the reverse sorted fnames by scores
	sortInd = [i[0] for i in sorted(enumerate(fnameScores),
	                                key = lambda x:x[1], reverse = True)]
	
	#print sortInd
	return [fnames[i] for i in sortInd[:10]]


def getTopTenProtsAllTargs(inDir, topProtDict, opt = 'avgMinDiff'):
	# Creates a directory for each 3mer target which holds 
	# the 10 most specific pwms/logos for that target using 
	# the optimization given by the opt parameter

	top10Dir = inDir + 'top10' + opt + '/'
	makeDir(top10Dir)
	top10Dict = {}

	for targ in sorted(topProtDict.keys()):
		print targ
		targPwmDir = inDir + targ + '/pwms/'
		targLogoDir = inDir + targ + '/logos/'
		handle = os.popen('ls ' + targPwmDir)
		nucMatDict = {}
		for fname in handle:
			fname = fname.strip()
			nucMatDict[fname] = pwmfile2matrix(targPwmDir + fname)
		
		# Get the list of top 10 pwm file names in decreasing order 
		top10 = getTopTenProts(targ, nucMatDict, opt)
		#print top10

		# Place these top 10 logos in a special directory
		top10TargDir = top10Dir + targ + '/'
		#print top10TargDir
		top10pwmDir = top10TargDir + 'pwms/'
		top10logoDir = top10TargDir + 'logos/'
		makeDir(top10TargDir)
		makeDir(top10pwmDir)
		makeDir(top10logoDir)


		for i, fname in enumerate(top10):
			new_fname = '_'.join([str(i).zfill(2), fname])
			logo_fname = fname.split('.')[0] + '.pdf'
			new_logo_fname = '_'.join([str(i).zfill(2), logo_fname])
			os.system('cp %s %s' %(targPwmDir + fname, top10pwmDir + new_fname))
			os.system('cp %s %s' %(targLogoDir + logo_fname, \
			          			   top10logoDir + new_logo_fname))

		# Record the best 10 logos into a dictionary
		top10Dict[targ] = ['_'.join(i.split('.')[0].split('_')[1:]) \
			for i in top10]

	# Write the dictionary to a nice spreadsheet
	fout = open(top10Dir + 'allTargBest10.txt', 'w')
	fout.write('\t'.join(["targ"] + [str(i + 1) for i in range(10)]) + '\n')
	for targ in top10Dict.keys():
		fout.write('\t'.join([targ] + top10Dict[targ]) + '\n')

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
	#for k in sorted(topProtDict.keys()):
	#	print k, len(topProtDict[k])

	# Create logos/pfms for each of the prots in topProtDict
	#predictTopProts(freqDict, topProtDict, outDir)

	# Find the top 10 most specific pwms for each 3mer target
	getTopTenProtsAllTargs(outDir, topProtDict, opt = 'minMinDiff') 
	

if __name__ == '__main__':
	main()