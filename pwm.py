# A module for creating pwm files.

import os
import re
import numpy as np
from scipy import stats

nucs = ['A', 'C', 'G', 'T']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
	      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Used by other modules

def infoEntr(f):
	# Returns the information entropy of the 
	# distribution f (base 2), using the convention
	# that plg(p) = 0 when p = 0.

	ent = 0.0
	for i in range(len(f)):
		if f[i] == 0:
			continue
		ent -= f[i] * np.log2(f[i])

	return ent

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
	
def getConsensus(nucmat):
	# Returns a list of consensus bases, one for 
	# each position, of a matrix of nucleotide 
	# frequencies for a binding site.  A base 
	# is considered a consensus base if it has
	# > 0.5 frequency.  Otherwise, if there is 
	# a pair that has combined freq of >0.75, this
	# is considered a consensus pair.
	# If no consensus or consensus pair can be 
	# found, the ouptut is 'X'
	
	nucs = ['A', 'C', 'G', 'T']
	consList = []
	for i in range(len(nucmat)):
		col = nucmat[i]

		# Is there a single consensus letter?
		maxFreq = np.max(col)
		if maxFreq > 0.5:
			consList.append(nucs[np.argmax(col)])
			continue

		# Is there a consensus pair?
		maxFreq1 = 0
		maxFreq1pos = -1
		maxFreq2 = 0
		maxFreq2pos = -1
		
		# Get the top 2 frequencies and positions
		for j in range(len(col)):
			freq = col[j]
			if freq > maxFreq1:
				maxFreq2 = maxFreq1
				maxFreq2pos = maxFreq1pos
				maxFreq1 = freq
				maxFreq1pos = j
			elif freq > maxFreq2:
				maxFreq2 = freq
				maxFreq2pos = j
		
		if maxFreq1 + maxFreq2 > 0.75:
			consList.append(nucs[maxFreq1pos] + nucs[maxFreq2pos])
		else:
			consList.append('X')

	return consList


def comparePCC(predPWM, expPWM):
	# This function asumes that the two PWMs 
	# are already aligned to one another and 
	# are of the same length.
	#
	# It scores the columns based on the information
	# content weighted Pearson correlation, adding 
	# the following for each column, i:
	# PCC(pred_i, exp_i) * IC(exp_i)/2
	colPcc = []
	colPcc_ic = []
	for i in range(len(predPWM)):
		corr = stats.pearsonr(predPWM[i,:], expPWM[i,:])
		ic = 2 - infoEntr(expPWM[i,:])
		colScore = corr[0] * ic/2
		colPcc.append(corr[0])
		colPcc_ic.append(colScore)
	
	#print colPcc, colPcc_ic	
	return colPcc, colPcc_ic

def pwmfile2matrix(pwmFile):
	# Converts a file of the transfac format to 
	# a matrix representaiton of the file
	fin = open(pwmFile, 'r')
	fin.readline()  # Skip the headers
	fin.readline()
	numLetters = len(fin.readline().strip().split()) - 1

	numList = []
	line = fin.readline()
	while line.strip() != 'XX':
		numList.append(line.strip().split()[1:-1])
		line = fin.readline()

	pwm = np.zeros((len(numList), numLetters), float)
	for i in range(len(numList)):
		pwm[i,:] = numList[i]

	return pwm

def makeLogo(infile, outfile, format = 'pdf', alpha = 'protein',
             composition = 'none', size = 'large', 
             colScheme = 'chemistry', xlab = None,
             ylab = None, annot = None, fineprint = None):
	# Creates a logo based on a pwm by calling weblogo.
	#
	# - infile is a path to a pwm file in the transfac 
	# file format.
	# - outfile is the path where the file will be saved
	# The remaining parameters correspond ot flags for weblogo.
	# See weblogo CLI for more info about them.
	#print "Making logo"
	opts = '-F %s -A %s -i 0 -s %s -c %s --composition %s' \
			%(format, alpha, size, colScheme, composition)
	if xlab != None:
		opts += ' -x %s' %xlab
	if annot != None:
		opts += ' --annotate %s' %annot
	if fineprint != None:
		opts += ' --fineprint %s' %fineprint
	if ylab != None:
		opts += ' --ylabel %s' %ylab
	#print opts
	#print infile
	#print outfile
	os.system('weblogo %s < %s > %s' %(opts, infile, outfile)) 
	pass



## Specific to this program ... not used by other modules.

def writePWM(dst, posCounts, npos, letters):
	# Output a PWM file based on a posCount
	# dictionary (pos, letter) -> count/freq

	fout = open(dst, 'w')

	# Normalize posCounts to a distribution
	for i in range(npos):
		rowTot = 0
		for j in letters:
			rowTot += posCounts[i,j]	
		if rowTot == 0:
			fout.close()
			os.remove(dst)
			return
		for j in letters:
			posCounts[i,j] = posCounts[i,j]/float(rowTot)

	# Write a dummy header
	fout.write('ID idNum\nBF species\n')
	fout.write('P0\t' + '\t'.join(letters) + '\n')
	for i in range(npos):
		total = 0
		outstr = str(i+1).zfill(2)
		for j in letters:
			outstr += '\t' + '%.4f' %(posCounts[i,j])
			total += posCounts[i,j]
		outstr += '\tX\n'
		fout.write(outstr)
	fout.write('XX\n\\\\\n')

def initPosCounts(npos, type):
	# Initialize the dictionay to all 0 counts
	nucs = ['A', 'C', 'G', 'T']
	aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
		      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	
	posCounts = {}  # (pos, letter) -> count
	for p in range(npos):
		if type == 'protein':
			for a in aminos:
				posCounts[(p,a)] = 0.0
		elif type == 'dna':
			for n in nucs:
				posCounts[(p,n)] = 0.0

	return posCounts

def parseCountFile(src, posCounts):
	# Parses the file of sequences and counts

	fin = open(src, 'r')
	for line in fin:
		sp_line = line.strip().split()
		seq, freq = sp_line[0], eval(sp_line[1])
		for i, letter in enumerate(seq):
			posCounts[i,letter] += freq

def makePwm(src, dst, npos, type):
	# src is a tab delimited txt file where the first column 
	# is a list of proteins or dnas (each of same length)
	# and the second is a list of counts or frequencies.
	# Outputs a pwm file in the transfac format to dst.
	#
	# npos is the number of positions per sequence
	# type is either 'protein' or 'dna'

	posCounts = initPosCounts(npos, type)
	letters = list(set(sorted([i[1] for i in posCounts.keys()])))
	#print letters

	# Fill the dictionary
	parseCountFile(src, posCounts)	

	# Output a PWM file in transfac format
	writePWM(dst, posCounts, npos, letters)

def makePWMDir(prefix):
	# Make pwms and logos for each of the sequence
	# files in the dirctory 'prefix'

	# Create the new directories
	print prefix
	try:
		os.mkdir(prefix + 'pwm/')
	except OSError:
		pass
	try:
		os.mkdir(prefix + 'logos/')
	except OSError:
		pass

	# For each seq file in the file
	handle = os.popen('ls ' + prefix)
	for fname in handle:
		fname = fname.strip()
		if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
			continue
		print "Processing %s" %(prefix + fname)
		targ = fname.split('.')[0]
		src = prefix + fname
		dst = prefix + 'pwm/' + fname
		makePwm(src, dst, 6, "protein")
		logoFname = prefix + 'logos/' + targ + '.pdf'
		makeLogo(dst, logoFname, xlab = targ)

def main():
	prefix = '../data/b1hData/newDatabase/6varpos/' + \
		'F2/high/protein_seq_cut10bc_0/'
	makePWMDir(prefix)

if __name__ == '__main__':
	main()