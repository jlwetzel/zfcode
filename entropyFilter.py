import os
import jellyfish
import numpy as np
from fixTables import normalizeFreq

def rmBadSeqs(fname, canonIndex):
	# Remove the sequences that have only one coding 
	# posibility and are more than edit distance 1 
	# away from all canonical sequences for this target 
	# that meet the entropy cutoff.
	
	# Get all sequences meeting the entropy threshold
	fin = open(fname, 'r')
	goodSeqs = [i.strip().split('\t')[0] for i in fin.readlines()\
	            if eval(i.strip().split('\t')[3]) != 1]
	fin.close()
	
	# Convert to canonical set of sequences
	canonSeqs = []
	for seq in goodSeqs:
		canonSeq = ''
		for i in canonIndex:
			canonSeq += seq[i]
		canonSeqs.append(canonSeq)
	canonSeqs = set(canonSeqs)
	#print canonSeqs

	# Write all the good sequences from fname into a temp 
	# file and then replace fname with the tmp file
	fin = open(fname, 'r')
	fout = open('tmp.txt', 'w')
	for line in fin:
		sp_line = line.strip().split('\t')
		numPoss = eval(sp_line[3])
		if numPoss > 1:
			fout.write(line)
		else:
			seq = sp_line[0]
			seq2 = ''
			for i in canonIndex:
				seq2 += seq[i]
			for canonSeq in canonSeqs:
				if jellyfish.hamming_distance(seq2, canonSeq) <= 1:
					#print "Here"
					fout.write(line)
					break
	fout.close()
	normalizeFreq('tmp.txt', 1)
	os.system('mv tmp.txt ' + fname)

def parseAllFile(prefix, rest, cutoff, canonIndex):
	# Create a set of new files with that have 
	# been filtered for entropy in a new directory.

	fin = open(prefix + rest, 'r')
	newdir = prefix + 'protein_cut10_entr' + \
	         str(cutoff).replace('.', '') + '/'
	try:
		os.mkdir(newdir)
	except OSError:
		pass

	lastTarg = 'AAA'
	fout = open(newdir + lastTarg + '.txt', 'w')
	print 'Processing '+ newdir + lastTarg + '.txt'
	i = 0
	j = 0
	for line in fin:
		sp_line = line.strip().split('\t')
		targ = sp_line[0]
		numObs = eval(sp_line[3])
		numPoss = eval(sp_line[4])
		entropy = eval(sp_line[6])

		if targ != lastTarg:
			fout.close()
			if i == 0:
				os.system('rm ' + newdir + 'AAA.txt')
				fout = open(newdir + targ + '.txt', 'w')
			else:
				print "%d unique sequences" %j
				rmBadSeqs(newdir + lastTarg + '.txt', canonIndex)
				j = 0
				fout = open(newdir + targ + '.txt', 'w')
				print 'Processing '+ newdir + targ + '.txt'
		lastTarg = targ
			
		if entropy >= cutoff or (numPoss == 1 and numObs == 1):
			fout.write('\t'.join(sp_line[1:]) + '\n')
			j += 1
		i += 1
	fout.close()
	rmBadSeqs(newdir + targ + '.txt', canonIndex)
	fin.close()

def main():

	npos = 5

	# For the 6 position data
	if npos == 6:
		cutoffs = {('F1', 'high'): 0.25, ('F1', 'low'): 0.25,
			   	   ('F2', 'high'): 0.25, ('F2', 'low'): 0.25,
			   	   ('F3', 'high'): 0.25, ('F3', 'low'): 0.25,}
		fings = ['F1', 'F2', 'F3']
		strins = ['low', 'high']
		tag = '6varpos'
		canonIndex = [0,2,3,5]

	# For the 5 position data
	elif npos == 5:
		cutoffs = {('F2', 'high'): 0.25, ('F2', 'low'): 0.25}
		fings = ['F2']
		strins = ['low', 'high']
		tag = '5varpos'
		canonIndex = [0,1,2,4]

	for f in fings:
		for s in strins:
			prefix = '../data/b1hData/newDatabase/' + tag \
				+ '/' + f + '/' + s + '/'
			rest = 'protein_seqs_JSD/all_cut10.txt'
			parseAllFile(prefix, rest, cutoffs[f,s], canonIndex)

if __name__ == '__main__':
	main()