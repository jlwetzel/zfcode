import os
import re
from jellyfish import hamming_distance
from fixTables import normalizeFreq
from gatherBindStats import getProtSet

def hasSupport(seq, keepSeqs, canInd):
	# Returns True if seq is at hamming 
	# distance of one from at least one 
	# other seq in keepSeqs when using 
	# only the canonical indices

	canSeq1 = ''
	for j in canInd:
		canSeq1 += seq[j]

	for s in keepSeqs:

		# Make sure this is not the exact same sequence
		if seq == s:
			continue

		canSeq2 = ''
		for j in canInd:
			canSeq2 += s[j]
		if hamming_distance(canSeq1, canSeq2) <= 1:
			return True

	return False


def filterEntropySupport(oldFile, newFile, cutoff, canInd):
	# Scans oldFile for lines that don't pass entropy cutoff.
	# These lines are removed, except in the case of proteins
	# with only one possible coding combination.  
	# The remaining lines will then be scanned once more to 
	# see if they have "support".  
	# Support means that there is at least one other sequence
	# that is hamming distance at most 1 from the sequence in 
	# the canonical positions -1,2,3,6.

	# Get list of lines/seqs we could potentially keep
	# (either pass entropy filter or only one way of coding)
	print "Processing %s" %oldFile
	fin = open(oldFile, 'r')
	keepLines = []
	for line in fin:
		sp_line = line.strip().split()
		entropy = eval(sp_line[5])
		numPoss = eval(sp_line[3])
		if entropy >= cutoff or numPoss == 1:
			keepLines.append(line)
	fin.close()

	# Only keep seqs from keepLines that are supported 
	# by at least one other seq in keepLines
	fout = open(newFile, 'w')
	keepSeqs = [i.strip().split()[0] for i in keepLines]
	for i, seq in enumerate(keepSeqs):
		if hasSupport(seq, keepSeqs, canInd):
			fout.write(keepLines[i])	
	fout.close()

def main():
	
	varpos = '5varpos'
	fings = ['F2', 'F3']
	strins = ['low', 'high']
	cutoff = 0.25
	for f in fings:
		for s in strins:
			oldDir = '../data/b1hData/newDatabase/' + \
				'/'.join([varpos, f, s, 'protein_seq']) + '/'

			# Make a new directory for the filtered proteins
			newDir = '/'.join(oldDir.split('/')[:-2]) + \
				'/protein_seq_'+ (str(cutoff)).replace('.', '') + '/'
			try:
				os.mkdir(newDir)
			except OSError:
				pass

			# Step through each protein file
			handle = os.popen('ls ' + oldDir)
			for fname in handle:
				fname = fname.strip()
				if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
					continue

				if varpos == '6varpos':
					canInd = [0,2,3,5]
				elif varpos == '5varpos':
					canInd = [0,1,2,4]
				filterEntropySupport(oldDir + fname, newDir + fname,
				                     cutoff, canInd)
				normalizeFreq(newDir + fname, 1)

if __name__ == '__main__':
	main()