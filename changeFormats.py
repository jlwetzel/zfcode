# A script for switching between ways of formatting the b1h data.

import os
import sys

nucs = ['A', 'C', 'G', 'T']

def targetFiles2CSV(dirPath):
	# Just puts a directory of nucleotide triplet binding files 
	# (NNN.txt) into one a single csv file where columns
	# are:
	# count n1 n2 n3 a0 a1 a2 a3 a5 a6
	#
	# Assumes files are named <NNN.txt> and in path dir
	# Does not remove any files.

	allFile = open(dirPath + "all.csv", 'w')
	allFile.write("count,n1,n2,n3,a0,a1,a2,a3,a5,a6\n")
	for n1 in nucs:
		for n2 in nucs:
			for n3 in nucs:
				triplet = n1+n2+n3
				try:
					tripletFile = open(dirPath + triplet+'_amino_acid_analysis.txt', 'r')
				except IOError:
					print "Invalid File Name: "+triplet+'_amino_acid_analysis.txt'
					pass
				for line in tripletFile:
					obs = ""
					prot = line.strip().split()[0]
					count = line.strip().split()[1]
					obs += count + ","
					for n in triplet:
						obs += n + ","
					for a in prot[:-1]:
						obs += a + ","
					obs += a[-1] + "\n"
					allFile.write(obs)
	allFile.close()

def csv2oneFile(dirPath):
	# Converts the csv file 'all.csv' inside dirPath to a 
	# plaintext file 'all.txt' with 3 columns (NNN, AAAAAA, Count).
	# Does not remove any files.
	csvFile = open(dirPath + 'all.csv', 'r')
	outfile = open(dirPath + 'all.txt', 'w')
	for line in csvFile:
		sp_line = line.strip().split(',')
		count = sp_line[0]
		base = sp_line[1] + sp_line[2] + sp_line[3]
		prot = sp_line[4] + sp_line[5] + sp_line[6] + sp_line[7] \
			+ sp_line[8] + sp_line[9]
		outfile.write(base + "\t" + prot + "\t" + count + "\n")
	csvFile.close()
	outfile.close()

def main():
	#path = sys.argv[1]
	targetFiles2CSV("../data/b1hData/oldData/F3/unfiltered/high/")
	csv2oneFile("../data/b1hData/oldData/F3/unfiltered/high/")


if __name__ == '__main__':
	main()