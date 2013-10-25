import os
import re
from jellyfish import hamming_distance
from fixTables import normalizeFreq
from gatherBindStats import getProtSet

def hasSupport(seq, keepSeqs, codedOneWay,
               supportCutoff, canInd):
	# Returns True if seq is at hamming 
	# distance of at most one from at least <supportCutoff> 
	# other (non-identical) seqs in keepSeqs when 
	# looking in canonical positions (-1,2,3,6)
	# and which has more than one possible coding variant

	if supportCutoff == 0:
		return True

	canSeq1 = ''
	for j in canInd:
		canSeq1 += seq[j]

	supportFound = 0
	for i, s in enumerate(keepSeqs):

		# Make sure this is not the exact same sequence
		if seq == s:
			continue

		# If find enough canonical prots within ham dist 1
		# that each have more than one possible coding variant
		# then return True
		canSeq2 = ''
		for j in canInd:
			canSeq2 += s[j]
		if hamming_distance(canSeq1, canSeq2) <= 1 \
			and not codedOneWay[i]:
			supportFound += 1
			if supportFound == supportCutoff:
				return True

	# Did not find enough support sequences
	return False


def filterEntropySupport(oldFile, newFile, freqCutoff, entCutoff, 
                         supportCutoff, canInd):
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
	codedOneWay = []
	for line in fin:
		sp_line = line.strip().split()
		#entropy = eval(sp_line[5])  #For my format
		#numPoss = eval(sp_line[3])  #For my format
		entropy = eval(sp_line[4])   #Anton's format
		#print entropy
		numPoss = eval(sp_line[2])   #Anton's format
		freq = eval(sp_line[1])      #Anton's format
		if freq >= freqCutoff and (entropy >= entCutoff or numPoss == 1):
			keepLines.append(line)
			if numPoss == 1:
				codedOneWay.append(True)
			else:
				codedOneWay.append(False)
	fin.close()

	# Only keep seqs from keepLines that are supported 
	# by at least one other seq in keepLines
	fout = open(newFile, 'w')
	keepSeqs = [i.strip().split()[0] for i in keepLines]
	for i, seq in enumerate(keepSeqs):
		if hasSupport(seq, keepSeqs, codedOneWay,
		              supportCutoff, canInd):
			fout.write(keepLines[i])	
	fout.close()

def main():
	
	varpos = '6varpos'
	fings = ['F1']#['F2', 'F3']
	strins = ['low', 'high']
	entCutoff = 0.25
	supportCutoff = 1
	freqCutoff = 0.0001
	oldProts = 'unfiltered2'

	for f in fings:
		for s in strins:
			oldDir = '../data/b1hData/antonProcessed/' + \
				'/'.join([f, s, oldProts]) + '/'
			
			# For Anton's data
			newDir = '/'.join(oldDir.split('/')[:-2]) + \
				'/' + 'filt_10e-4_'+ \
				(str(entCutoff)).replace('.', '') + \
				'_' + str(supportCutoff) + '_c/'
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
				                     freqCutoff,
				                     entCutoff, 
				                     supportCutoff, 
				                     canInd)
				normalizeFreq(newDir + fname, 1)

if __name__ == '__main__':
	main()