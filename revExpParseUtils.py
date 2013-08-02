import numpy as np

def get700Prots(fname):
	# Parses the file with all 700 proteins and
	# returns a tuple list of the for (number, F2 protein)

	print fname
	fin = open(fname, 'r')
	fin.readline()   # Skip the header line
	tList = []
	for line in fin:
		print line
		sp_line = line.strip().split()
		num = int(sp_line[0])
		prot = sp_line[2]
		tList.append((num,prot))
	fin.close()

	return tList

def normalizePWM(pwm):
	# Normalizes each inner array of 
	# a 2d numpy array
	for i in len(pwm):
		pwm[i,:] = pwm[i,:]/np.sum(pwm[i,:])
	return pwm

def parseMemeFile(fpath):
	# Returns a PWM of the width of the alignment in the
	# MEME.txt file.
	nucs = ['A', 'C', 'G', 'T']

	fin = open(fpath, 'r')
	# Find the alignments width
	line = ''
	for line in fin:
		print line
		if line.strip() == '':
			continue
		sp_line = line.strip().split()
		if sp_line[0] == 'MOTIF':
			alWid = eval(sp_line[4])
		if line.strip() == 'Motif 1 sites sorted by position p-value':
			break
	fin.readline()
	fin.readline()
	fin.readline()

	# Build the count per-position count matrix
	pwm = np.zeros((alWid, len(nucs)), float)

	for line in fin:
		sp_line = line.strip().split()
		count = np.log2(eval(sp_line[0].split('_')[1]))
		nucseq = sp_line[4]
		for i, n in enumerate(nucseq):
			pwm[i,nucs.index(n)] += count

	return normalizePWM(pwm)

def makeallpwms700s(listfile): 
	# Makes pwms for all of the proteins that end 
	# in a 700 (Think this is all F2s?)
	bcKey = {1: 'TC', 2: 'AA', 3: 'GG'}
	fin = open(listfile, 'r')
	for line in fin:
		sp_line = line.strip().split()
		bc = bcKey[eval(sp_line[4][-1])] + '_' + sp_line[2]
		dset = sp_line[0].split('-')[0]
		fpath = '../data/revExp/MN28-37/%s/%s/meme.txt' \
			%(dset, bc)
		pwm = parseMemeFile(fpath)
		print pwm


def main():
	fname = '../data/revExp/revExpBarcodes/all700Entries.txt'
	makeallpwms700s(fname)

if __name__ == '__main__':
	main()