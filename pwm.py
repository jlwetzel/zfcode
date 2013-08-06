# A module for creating pwm files.

import os
import re
import numpy as np
from scipy import stats

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

def comparePWMs(predPWM, expPWM, threshold = 0.25):
	# This function asumes that the two PWMs 
	# are already aligned to one another and 
	# are of the same length.
	#
	# It scores the columns based on the information
	# content weighted Pearson correlation, adding 
	# the following for each column, i:
	# PCC(pred_i, exp_i) * IC(exp_i)/2
	score = 0
	numCorrect = 0
	numCorrectIC = 0
	for i in range(len(predPWM)):
		corr = stats.pearsonr(predPWM[i,:], expPWM[i,:])
		ic = 2 - infoEntr(expPWM[i,:])
		colScore = corr[0] * ic/2
		if corr[0] >= threshold:
			numCorrect += 1
		if colScore >= threshold:
			numCorrectIC += 1
		score += colScore
		
	return score, numCorrect, numCorrectIC, len(predPWM)

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
             annot = None, fineprint = None):
	# Creates a logo based on a pwm by calling weblogo.
	#
	# - infile is a path to a pwm file in the transfac 
	# file format.
	# - outfile is the path where the file will be saved
	# The remaining parameters correspond ot flags for weblogo.
	# See weblogo CLI for more info about them.
	opts = '-F %s -A %s -i 0 -s %s -c %s --composition %s' \
			%(format, alpha, size, colScheme, composition)
	if xlab != None:
		opts += ' -x %s' %xlab
	if annot != None:
		opts += ' --annotate %s' %annot
	if fineprint != None:
		opts +  ' --fineprint %s' %fineprint
	#print opts
	#print infile
	#print outfile
	os.system('weblogo %s < %s > %s' %(opts, infile, outfile)) 
	pass

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
	# Make pwms and logos for each of the files sequence
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
		'F2/low/protein_seq_cut3bc_025/'
	makePWMDir(prefix)

if __name__ == '__main__':
	main()