# A binding specificity predictor based on log-likelihood ratios
# for the important contacts.  
#
# Requires first that a log-likelihood txt table had been created 
# by functions from mlUtils.py and ./figurecode/llAnalysis.R modules

import numpy as np
import os
from pwm import makeLogo, pwmfile2matrix, comparePCC, makeNucMatFile, getConsensus
from revExpParseUtils import getTargDict
import re
from itertools import combinations_with_replacement

nucs = ['A', 'C', 'G', 'T']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
	      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Maps base or helix positions (0 corresponds to -1 for helix)
# to indicies for canonical helices
baseposMap = {1:0, 2:1, 3:2}
aminoposMap = {0:0, 2:1, 3:2, 6:3}

class LLPredictor():

	def __init__(self, table):
		# table is a path to the table file created by llAnalysis.R
		fin = open(table, 'r')
		headerLine = fin.readline()
		self.contacts = headerLine.strip().split()[1:]
		self.contacts = [(i[1], i[3]) for i in self.contacts]
		self.llmats = self.makeLLMatsDict(fin)
		fin.close()
		

	def makeLLMatsDict(self, fin):
		# Create a dictionary of log-likelihood matrices
		llmats = {}
		for line in fin:
			sp_line = line.strip().split()
			c = (eval(sp_line[0][1]), eval(sp_line[0][3]))
			b = sp_line[1]
			a = sp_line[2]
			aInd = aminos.index(a)
			bInd = nucs.index(b)
			score = eval(sp_line[3])
			if llmats.has_key(c):
				llmats[c][bInd,aInd] = score
			else:
				llmats[c] = np.zeros((len(nucs), len(aminos)),
				                      dtype = 'float')
				llmats[c][bInd,aInd] = score

		#llmats = self.centerMatColumns(llmats)
		return llmats

	def centerMatColumns(self, llmats):
		# Mean center columns of the matrix
		for k in llmats.keys():
			for j in range(len(aminos)):
				meanVal = np.mean(llmats[k][:,j])
				llmats[k][:,j] = llmats[k][:,j] - meanVal
		return llmats

	def predictCanonZF(self, prot):
		# Predicts the 3bp nucleotide specificty of prot, 
		# which is are the canonical (-1,2,3,6) positions
		# of the alpha helix of a ZF domain

		temp = 2   # temperature parameter for boltzmann distr
		nucMat = np.zeros(shape = (3,4), dtype = float)
		
		# Make a dictionary of lists of amino positions positions
		# to predict bases from ... bPos -> [aPos1, aPos2 ..., aPos_n]
		toPredict = {}
		for k in sorted(self.llmats.keys()):
			bpos = k[0]
			apos = k[1]
			if toPredict.has_key(bpos):
				toPredict[bpos].append(apos)
			else:
				toPredict[bpos]= [apos]
		#print toPredict

		# Get the base frequency vector for each 
		# position to be predicted
		for bpos in sorted(toPredict.keys()):
			for apos in toPredict[bpos]:
				a = prot[aminoposMap[apos]]
				#print a
				os_bpos = baseposMap[bpos]
				#print bpos, apos, self.llmats[(bpos, apos)]
				nucMat[os_bpos] += \
					self.llmats[(bpos, apos)][:,aminos.index(a)]
			
			# Average if there were multiple contacts for the base
			nucMat[os_bpos] = nucMat[os_bpos]/float(len(toPredict[bpos]))

			# Convert back to get prob ratio and normalize
			nucMat[os_bpos] = np.exp(temp*nucMat[os_bpos])
			nucMat[os_bpos] = nucMat[os_bpos]/np.sum(nucMat[os_bpos])

		return nucMat

	def predictCanonZFArray(self, canonZFs):
		# Returns a 2d numpy array of the predicted pwm for 
		# a list of arrayed ZFs.  Assumes that the each ZF
		# in canonZFs is simply a length 4 string corresponding
		# to positions -1,2,3,6 of the alpha helix.

		numZFs = len(canonZFs)
		pwm = np.zeros(shape = (numZFs*3, 4))
		
		for i in range(numZFs):
			nmat = self.predictCanonZF(canonZFs[i])
			for j in range(len(nmat)):
				pwm[i*len(nmat) + j,:] = nmat[j,:]
		
		return pwm

def testLLPredictor(table):
	pred = LLPredictor(table)
	#for k in sorted(pred.llmats.keys()):
	#	print k
	#	print pred.llmats[k]
	x = pred.predictCanonZFArray(["RDLR", "REHR"])
	print x
	label = "RDLR-REHR"
	makeNucMatFile('../data/scratch/', label, x)
	logoIn = '../data/scratch/' + label + '.txt'
	logoOut = '../data/scratch/' + label + '.pdf'
	makeLogo(logoIn, logoOut, alpha = 'dna', 
			         colScheme = 'classic',
			         #annot = "'5,M,3'",
			         xlab = label)

def predictMarcusPWMs(predictor, outputDir, finger, strin, 
                      filt, pred = 'llratio'):
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
	fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
	           %('num', 'targ','prot','canonprot','c1pcc','c2pcc',
	             'c3pcc', 'c1pcc.ic', 'c2pcc.ic', 'c3pcc.ic', 'p.c1cons',
	             'p.c2cons', 'p.c3cons', 'e.c1cons', 'e.c2cons', 'e.c3cons',
	             'pred.filt', 'finger', 'strin'))

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
		goal = sp_fname[1]
		prot = sp_fname[2].split('.')[0]
		canonProt = prot[0] + prot[2] + prot[3] + prot[6]
		label = '_'.join([str(protNum), goal, prot, strin])
		
		# Make the prediciton and write out to the correct files.
		nmat = predictor.predictCanonZF(canonProt)
		makeNucMatFile(pwmdir, label, nmat)
		logoIn = pwmdir + label + '.txt'
		logoOut = logodir + label + '.pdf'
		makeLogo(logoIn, logoOut,
		         alpha = 'dna', colScheme = 'classic',
		         annot = "'5,M,3'",
		         xlab = '_'.join([goal,prot]))

		# Compare this pwm to the reverse experiment
		expMat = pwmfile2matrix(expDir + fname)
		colPcc, colPcc_ic = comparePCC(nmat, expMat)
		predCons = getConsensus(nmat)
		expCons = getConsensus(expMat)

		# Write the comparison results to file
		fout.write("%s\t%s\t%s\t%s\t" %(protNum, goal, prot, canonProt))
		fout.write("%.4f\t"*6 %(colPcc[0], colPcc[1], colPcc[2], \
		           				colPcc_ic[0], colPcc_ic[1], colPcc_ic[2]))
		fout.write("%s\t"*6 %(predCons[0], predCons[1], predCons[2],
		           			  expCons[0], expCons[1], expCons[2]))
		fout.write("%s\t%s\t%s\n" %(pred+'.'+filt, finger, strin))

	fout.close()

def makeDir(path):
	# Try to make a path and pass if it can't be created
	try:
		os.mkdir(path)
	except OSError:
		pass

def main():

	#testLLPredictor('../figures/llAnalysis/F2/low/cut10bc_0_5/llRatios_c.txt')
	
	testFing = 'F2'
	testStrin = 'low'
	fings = ['F2']#['F3', 'F2']
	strins = ['low']
	filts = ['cut10bc_0_5']#, 'cut3bc_0_5', 'cut3bc_025', 'cut10bc_025']
	filtsLabs = ['c10_0_5']#, 'c3_0_5', 'c3_025', 'c10_025']
	for f in fings:
		for s in strins:
			for i, filt in enumerate(filts):
				inDir = '../figures/llAnalysis/' \
					+ f + '/' + s + '/' + filt + '/'
				outDir = '../data/llRatioCA/' + f + '/' + s + \
					'/' + filt + '/'
				predictor = LLPredictor(inDir + 'llRatios_ca.txt')
				predictMarcusPWMs(predictor, outDir, testFing, testStrin,
				                  filtsLabs[i])
				

if __name__ == '__main__':
	main()