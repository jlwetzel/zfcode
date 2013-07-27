import os
import re
import numpy as np

def normalizeFreq(fname, coln):
	# Normalizes the numbers is column coln to a 
	# proper distribution.  column numbers start
	# at 0.
	# This is a fast method that assumes the whole 
	# table is large enough to fit in memory!
	fin = open(fname, 'r')
	fout = open('tmp2.txt', 'w')
	lines = [i.strip().split() for i in fin]
	fin.close()

	freqs = np.array([eval(i[coln]) for i in lines])
	freqs = list(freqs/float(np.sum(freqs)))
	for i, line in enumerate(lines):
		line = '\t'.join(line[:coln] + [str(freqs[i])] + \
		                 line[coln+1:]) + '\n'
		fout.write(line)
	fout.close()
	os.system('mv ' + 'tmp2.txt ' + fname)

def remLowCount(fname, coln, cutoff):
	# Removes any line if the file if the 
	# value in coln for that line <= cutoff
	fin = open(fname, 'r')
	fout = open('tmp2.txt', 'w')
	for line in fin:
		sp_line = line.strip().split()
		if eval(sp_line[coln]) > cutoff:
			fout.write(line)
	fin.close()
	fout.close()
	os.system('mv ' + 'tmp2.txt ' + fname)

def remLowCountDir(dirname, coln, cutoff):
	# Removes any line if the file if the 
	# value in coln for that line <= cutoff
	# for every binding file in the directory
	handle = os.popen('ls ' + dirname)
	for fname in handle:
		fname = fname.strip()
		if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
			continue
		remLowCount(dirname + fname, coln, cutoff)

def normalizeDir(dirname, coln):
	# Normalizes the numbers is column coln to a 
	# proper distribution.  column numbers start
	# at 0.  Does so for every binding file in 
	# the directory.
	handle = os.popen('ls ' + dirname)
	for fname in handle:
		fname = fname.strip()
		if re.match(r'[ACGT]{3}(.)*.txt', fname) == None:
			continue
		normalizeFreq(dirname + fname, coln)


def main():
	fing = 'F2'
	strin = 'high'
	protDir = 'threshold025'
	dirpath = '../data/b1hData/oldData/' + \
		'/'.join([fing, protDir, strin]) + '/'
	print dirpath
	#remLowCountDir(dirpath, 1, 10)
	#normalizeDir(dirpath, 1)

if __name__ == '__main__':
	main()