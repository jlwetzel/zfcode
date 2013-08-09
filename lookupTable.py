import os
import re
import numpy as np
from pwm import makeLogo, pwmfile2matrix, comparePWMs, makeNucMatFile
from fixTables import normalizeFreq
#from gatherBindStats import getProtDict


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
                     canonical, canInd, neighbors):
	# Appends a tuple to targList if a nearest
	# neighbor to the protein is found when 
	# allowing the 1, 3, or 6 positions of 
	# the alpha helix to vary.
	# Assumes files named like AAA.txt	

	targ = fname.split('/')[-1].split('.')[0]
	fin = open(fname, 'r')
	
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
			
		# Add to the frquency if this binder is a neighbor
		if binder in neighbors:
			totFreq += freq

	# Append to the targList
	#print targ, totFreq
	targList.append([targ, totFreq])


def get3merList(dirpath, varpos, protein, canonical = False,
                useNN = False, skipExact = False):
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
			swapInd = [0,2,3]
		else:
			swapInd = [0,3,5]		
		for i in swapInd:
			for a in aminos:
				if a != protein[i]:
					neighbors.append(protein[:i]+a+\
					                 protein[i+1:])
		#print "Neighbors:"
		#print neighbors

		handle = os.popen('ls ' + dirpath, 'r')
		for fname in handle:
			fname = fname.strip()
			if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
				continue
			updateTargListNN(dirpath + fname, targList, 
		                     protein, canonical, canInd, 
		                     neighbors)


	# Normalize the frequencies across the bound 3mers to 1
	totFreq = 0.0
	for [targ, freq] in targList:
		totFreq += freq
	for i in targList:
		i[1] = i[1]/float(totFreq)

	return [(i[0], i[1]) for i in targList]

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

def lookupMarcusPWMs(inDir, outputDir, finger, strin, 
                     filt, pred, 
                     useNN = False, skipExact = False):
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
	
	# Right now only working with 6 variable position 
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
		
		# Make the list of targets bound and normalized frequencies.
		targList = get3merList(inDir, 6, canonProt, canonical,
		                       useNN, skipExact)
		# Do something else if the target list is empty
		if targList == []:
			continue
			# Apply nearest neighbor strategy here if targList 
			# is empty instead of just ignoring?
		
		# Convert the 3mer target list to a frequency matrix,
		# write to file, and create a logo
		nucMat = targListToFreqMat(targList)
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

	# Don't use nearest neighbors
	fings = ['F2']
	strins = ['low']
	filts = [ 'cut10bc_0', 'cut3bc_025', 'cut10bc_025']
	filtsLabs = ['c3_025', 'c10_025', 'c10_0']
	for f in fings:
		for s in strins:
			for i, filt in enumerate(filts):
				inDir = '../data/b1hData/newDatabase/6varpos/' \
					+ f + '/' + s + '/' + 'protein_seq_' + filt + '/'
				outDir = '../data/lookupTable/'	+ f + '/' + s + \
					'/' + filt + '/'
				
				lookupMarcusPWMs(inDir, outDir, f, s, filtsLabs[i],
				                 'look', useNN = False, skipExact = False)
	"""
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

	# Use nearest neighbors and skip all exact matches
	fings = ['F2']
	strins = ['low']
	filts = ['cut10bc_0', 'cut3bc_025', 'cut10bc_025']
	filtsLabs = ['c3_025', 'c10_025', 'c10_0']
	for f in fings:
		for s in strins:
			for i, filt in enumerate(filts):
				inDir = '../data/b1hData/newDatabase/6varpos/' \
					+ f + '/' + s + '/' + 'protein_seq_' + filt + '/'
				outDir = '../data/lookupTableNNonly/' + f + '/' + s + \
					'/' + filt + '/'
				lookupMarcusPWMs(inDir, outDir, f, s, filtsLabs[i],
				                 'look.nnOnly', useNN = True, skipExact = True)


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