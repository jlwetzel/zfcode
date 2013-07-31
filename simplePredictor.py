# Creates a simple predictor where we create for each possible 
# contact pair and nuc/amino pairing a score based on frequency of 
# observed interations.

import numpy as np
import os
from pwm import makeLogo

nucs = ['A', 'C', 'G', 'T']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
	      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


class SimplePredictor():
	def __init__(self, protFile, numNucs, numAminos):
		self.protFile = protFile
		self.numNucs = numNucs
		self.numAminos = numAminos
		self.pairMap = self.makePosPairMap()
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


def main():
	numAminos = 6
	proteinDir = '../data/b1hData/newDatabase/6varpos/' +\
		'F2/low/protein_seq_cut3bc_025/'
	outputDir = '../data/simplePredictor/cut3bc_025/'

	predictor = SimplePredictor(proteinDir + 'all.txt', 
	                            3, numAminos)
	#for i in sorted(predictor.pairMats.keys()):
	#	print predictor.pairMap[i]
	#	print predictor.pairMats[i]
	predictor.writePosSpecPWMs(outputDir)
	predictor.makePosSpecLogos(outputDir)

if __name__ == '__main__':
	main()

