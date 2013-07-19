# A script for switching between ways of formatting the b1h data.

import os
import sys

nucs = ['A', 'C', 'G', 'T']

def targetFiles2singleCSV(dirPath, numVarPos):
	# Just puts a directory of nucleotide triplet binding files 
	# (NNN.txt) into one a single csv file where columns
	# are:
	# freq n1 n2 n3 a0 a1 a2 a3 a5 a6 obsCode possCode jsd
	#
	# Assumes files are named <NNN.txt> and in path dir
	# Does not remove any files.

	allFile = open(dirPath + "all.csv", 'w')
	
	if numVarPos == 6:
		allFile.write('freq,n1,n2,n3,a0,a1,a2,a3,a4,a6,' + \
	    	          'obsCode,possCode,jsd,entropy\n')
	elif numVarPos == 5:
		allFile.write('freq,n1,n2,n3,a0,a2,a3,a4,a6,' + \
	    	          'obsCode,possCode,jsd,entropy\n')
	protSet = set()
	j = 0
	for n1 in nucs:
		for n2 in nucs:
			for n3 in nucs:
				triplet = n1+n2+n3
				try:
					tripletFile = open(dirPath + triplet+\
					                   '_protein_seqs_JSD.txt', 'r')
				except IOError:
					print "Invalid File Name: "+triplet+\
						'_protein_seqs_JSD.txt'
					continue
				i = 0
				for line in tripletFile:
					sp_line = line.strip().split('\t')
					obs = ""
					prot, freq = sp_line[0], sp_line[1]
					obsCode, possCode, jsd, entropy = sp_line[2], \
						sp_line[3], sp_line[4], sp_line[5]
					protSet.add(prot)
					
					obs += freq + ','
					for n in triplet:
						obs += n + ','
					for a in prot:
						obs += a + ','
					obs += ','.join([obsCode, possCode, jsd, entropy])
					obs += '\n'
					allFile.write(obs)
					i += 1
				print "Num proteins in %s_protein_seqs_JSD.txt:  %d" \
					%(dirPath + triplet, i)
				j += i
	print "Total num protein/dna observations in all files:  %d" %j
	print "Total num unique proteins in all files:  %d" %len(protSet)
	print
	allFile.close()

def csv2txtFile(dirPath, numVarPos):
	# Converts the csv file 'all.csv' inside dirPath to a 
	# plaintext file 'all.txt' with 3 columns (NNN, AAAAAA, Count).
	# Does not remove any files.
	csvFile = open(dirPath + 'all.csv', 'r')
	outfile = open(dirPath + 'all.txt', 'w')
	i = 0
	header = csvFile.readline()
	for line in csvFile:
		sp_line = line.strip().split(',')
		freq = sp_line[0]
		base = ''.join(sp_line[1:4])
		
		if numVarPos == 6:
			prot = ''.join(sp_line[4:10])
			obsCode = sp_line[10]
			possCode = sp_line[11]
			jsd = sp_line[12]
			entropy = sp_line[13]
		elif numVarPos == 5:
			prot = ''.join(sp_line[4:9])
			obsCode = sp_line[9]
			possCode = sp_line[10]
			jsd = sp_line[11]
			entropy = sp_line[12]

		outStr = '\t'.join(base, prot, freq, obsCode, possCode,
		                   jsd, entropy) + '\n'
		outfile.write(outStr)
		i += 1
	print "Num obs processed in csv to txt conversion:  %d" %i
	csvFile.close()
	outfile.close()

def main():
	pass
	#path = sys.argv[1]
	#targetFiles2CSV("../data/b1hData/oldData/F3/unfiltered/high/")
	#csv2oneFile("../data/b1hData/oldData/F3/unfiltered/high/")


if __name__ == '__main__':
	main()