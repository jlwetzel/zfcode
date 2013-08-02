# Creates a simple predictor where we create for each possible 
# contact pair and nuc/amino pairing a score based on frequency of 
# observed interations.

import numpy as np
import os
from pwm import makeLogo
#from revExpParseUtils import get700Prots

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

def makeNucMatFile(path, prot, protLabel, nucMat, a = 'dna'):
	# Write a transfac style frequency matrix to file

	fout = open(path + protLabel + '.txt', 'w')
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

def main():
	numAminos = 6
	proteinDir = '../data/b1hData/newDatabase/6varpos/' +\
		'F2/low/protein_seq_cut3bc_025/'
	outputDir = '../data/simplePredictor/cut3bc_025/'

	predictor = SimplePredictor(proteinDir + 'all.txt', 
	                            3, numAminos)
	
	# ONly need these lines once.
	#predictor.writePosSpecPWMs(outputDir)
	#predictor.makePosSpecLogos(outputDir)
	
	prots = get700Prots('../data/revExp/revExpBarcodes/' + \
	                    'revExper_GAG_700s.txt')
	print prots

	for p in prots:
		protNum = str(p[0])
		prot = ''.join([p[1][0],p[1][2],p[1][3],p[1][6]])
		protLabel = protNum + '-' + prot
		#print p, prot, protNum, protLabel

		nmat = predictor.predictCanon(prot)
		predicitonDir = outputDir+'predictions/'
		makeNucMatFile(predicitonDir + 'pwms/', prot,
		               protLabel, nmat)
		makeLogo(predicitonDir+ 'pwms/' + protLabel+'.txt', 
		         predicitonDir+ 'logos/' + protLabel+'.pdf',
		         alpha = 'dna', colScheme = 'classic',
		         annot = "'5,M,3'", xlab = prot)

if __name__ == '__main__':
	main()

