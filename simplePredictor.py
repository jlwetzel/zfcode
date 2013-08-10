# Creates a simple predictor where we create for each possible 
# contact pair and nuc/amino pairing a score based on frequency of 
# observed interations.

import numpy as np
import os
from pwm import makeLogo, pwmfile2matrix, comparePWMs, makeNucMatFile
from revExpParseUtils import getTargDict
import re

nucs = ['A', 'C', 'G', 'T']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
	      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


class SimplePredictor():
	def __init__(self, protFile, numNucs, numAminos):
		self.protFile = protFile
		self.numNucs = numNucs
		self.numAminos = numAminos
		self.pairMap = self.makePosPairMap()
		self.revPairMap = self.makeRevPairMap()
		self.numObs, self.freqs = self.makePositionSpecFreqs()
		self.notObserved = {}
		self.pairMats = self.makePairMats()
		
		
	def makePosPairMap(self):
		# Maps a dictionary x value to a nuc-amino positon pairing
		# E.g. posPairMap[0] = (1, -1)
		posPairMap = {}

		realNucPos = [1,2,3]
		if self.numAminos == 6:
			realAminoPos = [-1,1,2,3,4,6]
		elif self.numAminos == 5:
			realAminoPos = [-1,2,3,4,6] 
		for i, nval in enumerate(realNucPos):
			for j, aval in enumerate(realAminoPos):
				posPairMap[i*self.numAminos + j] = (nval, aval)

		return posPairMap

	def makeRevPairMap(self):
		# Make the map from actual position pairs to 
		# x values.
		vals = [i[1] for i in self.pairMap.items()]
		keys = [i[0] for i in self.pairMap.items()]
		revPairMap = {}
		for i, val in enumerate(vals):
			revPairMap[val] = keys[i]
		return revPairMap

	def makePositionSpecFreqs(self):
		# Returns a dictionary of position specific pairwise
		# frequencies for all the proteins in the 

		# Initialize the  frequency dictionary
		freqDict = {}
		for x in range(self.numNucs * self.numAminos):
			for n in nucs:
				for a in aminos:
					freqDict[x,n,a] = 0.0

		# Update the frequency dictionary per line in the all.txt file
		fin = open(self.protFile, 'r')
		numObs = 0
		for line in fin:
			sp_line = line.strip().split()
			nucSeq, aminoSeq = sp_line[0], sp_line[1]
			freq = eval(sp_line[2])
			for i, n in enumerate(nucSeq):
				for j, a in enumerate(aminoSeq):
					freqDict[i*self.numAminos + j, n, a] += freq
			numObs += 1
		fin.close()

		return numObs, freqDict

	def makePairMats(self):
		# Makes a set of 2d numpy arrays, one per pair
		# of nucleotide amino acid positions.
		# index 1 is the amino, index 2 is the nuc.
		# (Both are in alphabetical order.)
		# Frequency values are normalized per amino acid
		# for each matrix.

		pairMats = {}
		for p in self.pairMap.keys():
			pairMats[p] = np.zeros(shape = (len(aminos), len(nucs)))
			
			for i, a in enumerate(aminos):
				for j, n in enumerate(nucs):
					pairMats[p][i,j] = self.freqs[p,n,a]

			for p in pairMats.keys():
				for i in range(len(aminos)):

					# Deal with the case where a particular amino
					# acid is not observed at all for this 
					# position pair.
					if np.sum(pairMats[p][i,:]) == 0:
						
						# Make a not that this was not obseved
						if self.notObserved.has_key(p):
							self.notObserved[p].append(aminos[i])
						else:
							self.notObserved[p] = [aminos[i]]

						# Give 0 info content to this amino acid
						pairMats[p][i,:] = \
							np.array([0.25, 0.25, 0.25, 0.25])

					pairMats[p][i,:] /= np.sum(pairMats[p][i,:])

		return pairMats

	def writePosSpecPWMs(self, outputDir):
		# Writes the position-specific frequency matrices 
		# to the posSpecPWM directory in transfac format

		try:
			os.mkdir(outputDir + 'posSpecPWMs/')
		except OSError:
			pass

		for p in self.pairMats.keys():
			npos, apos = self.pairMap[p]
			pwmfname = outputDir + 'posSpecPWMs/' + \
				'posPair_b%da%d.txt' %(npos, apos)

			fout = open(pwmfname, 'w')
			# Write a dummy header
			fout.write('ID idNum\nBF species\n')
			fout.write('P0\t' + '\t'.join(nucs) + '\n')
			for i in range(len(self.pairMats[p])):
				outstr = str(i+1).zfill(2)
				for j in range(len(self.pairMats[p][i,:])):
					outstr += '\t%.4f' %self.pairMats[p][i,j]
				outstr += '\tX\n'
				fout.write(outstr)
			fout.write('XX\n\\\\\n')
			fout.close()

	def makePosSpecLogos(self, outputDir):
		# Creates the positions-specific logos from the 
		# transface matrices

		try:
			os.mkdir(outputDir + 'posSpecLogos/')
		except OSError:
			pass

		for p in self.pairMats.keys():
			npos, apos = self.pairMap[p]
			posStr = 'b%da%d' %(npos, apos)

			annotList = "'A,C,D,E,F,G,I,H,K,L,M,N,P,Q,R,S,T,V,W,Y'"
			pwmfname = outputDir + 'posSpecPWMs/' + \
				'posPair_%s.txt' %(posStr)
			logofname = outputDir + 'posSpecLogos/' + \
				'posPair_%s.pdf' %(posStr)
			makeLogo(pwmfname, logofname, alpha = 'dna',
			         colScheme = 'classic', annot = annotList,
			         xlab = posStr, fineprint = '')


	def predictCanon(self, prot):
		# Returns a 2d numpy array for the predicted logo
		# of the given protein.
		# Prot is assumed to be given as a length 4 
		# string representing the -1,2,3,6 positions 
		# of the alpha helix of the ZF domain.

		nucMat = np.zeros(shape = (3,4))
		
		nucMat[0,:] = \
			self.pairMats[self.revPairMap[(1,6)]]\
				[aminos.index(prot[3]),:]
		nucMat[1,:] = \
			self.pairMats[self.revPairMap[(2,3)]]\
				[aminos.index(prot[2]),:]
		nucMat[2,:] = \
			(self.pairMats[self.revPairMap[(3,-1)]]\
			 	[aminos.index(prot[0]),:]+\
			self.pairMats[self.revPairMap[(3,2)]]\
			 	[aminos.index(prot[1]),:])/2
		
		return nucMat

	def predictCanonArray(self, canonZFs):
		# Returns a 2d numpy array of the predicted pwm for 
		# a list of arrayed ZFs.  Assumes that the each ZF
		# in canonZFs is simply a length 4 string corresponding
		# to positions -1,2,3,6 of the alpha helix.

		numZFs = len(canonZFs)
		pwm = np.zeros(shape = (numZFs*3, 4))
		
		for i in range(numZFs):
			nmat = self.predictCanon(canonZFs[i])
			for j in range(len(nmat)):
				pwm[i*len(nmat) + j,:] = nmat[j,:]
		
		return pwm

def makeNucMatFile(path, label, nucMat):
	# Write a transfac style frequency matrix to file

	fout = open(path + label + '.txt', 'w')
	# Write a dummy header
	fout.write('ID idNum\nBF species\n')
	fout.write('P0\t' + '\t'.join(nucs) + '\n')
	for i in range(len(nucMat)):
		outstr = str(i+1).zfill(2)
		for j in range(len(nucMat[i,:])):
			outstr += '\t%.4f' %nucMat[i,j]
		outstr += '\tX\n'
		fout.write(outstr)
	fout.write('XX\n\\\\\n')
	fout.close()

def outputPosSpecPWMs(predictor, outputDir):
	# Output the position-pair logos
	predictor.writePosSpecPWMs(outputDir)
	predictor.makePosSpecLogos(outputDir)

def makeDir(path):
	# Try to make a path and pass if it can't be created
	try:
		os.mkdir(path)
	except OSError:
		pass

def predictMarcusPWMs(predictor, outputDir, finger, strin, 
                      filt, pred = 'simp'):
	# Make predictions for the F2 reverse experiments.

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

	# Try to do a simple prediction for each of the 
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
		targ = sp_fname[1]
		prot = sp_fname[2].split('.')[0]
		canonProt = prot[0] + prot[2] + prot[3] + prot[6]
		label = '_'.join([str(protNum), targ, prot, strin])
		
		# Make the prediciton and write out to the correct files.
		nmat = predictor.predictCanon(canonProt)
		makeNucMatFile(pwmdir, label, nmat)
		logoIn = pwmdir + label + '.txt'
		logoOut = logodir + label + '.pdf'
		makeLogo(logoIn, logoOut,
		         alpha = 'dna', colScheme = 'classic',
		         annot = "'5,M,3'",
		         xlab = '_'.join([targ,prot]))

		# Compare this pwm to the reverse experiment
		score, colcor, colcorIC, totCol = comparePWMs(nmat, 
		                pwmfile2matrix(expDir + fname))
		
		# Write the comparison results to file
		fout.write("%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%d\t%s\t%s\t%s\n" \
		           %(protNum, targ, prot, canonProt, score, colcor, \
		           colcorIC, totCol, pred+'.'+filt, finger, strin))			
	fout.close()

def main():
	numAminos = 6

	fings = ['F2']
	strins = ['low']
	filts = ['cut10bc_0_5', 'cut3bc_0_5'] #[ 'cut10bc_0', 'cut3bc_025', 'cut10bc_025']
	filtsLabs = ['c10_0_5', 'c3_0_5'] #['c3_025', 'c10_025', 'c10_0']
	for f in fings:
		for s in strins:
			for i, filt in enumerate(filts):
				inDir = '../data/b1hData/newDatabase/6varpos/' \
					+ f + '/' + s + '/' + 'protein_seq_' + filt + '/'
				outDir = '../data/SimplePredictor/'	+ f + '/' + s + \
					'/' + filt + '/'
				predictor = SimplePredictor(inDir + 'all.txt', 
	                            			3, numAminos)
				predictMarcusPWMs(predictor, outDir, f, s, filtsLabs[i])

	
	# Output position specific matrices/logos
	# Only need this line once.
	#outputPosSpecPWMs(predictor, outputDir)
	
	

	

if __name__ == '__main__':
	main()

