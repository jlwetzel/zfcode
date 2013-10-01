# A script for switching between ways of formatting the b1h data.

import os
import sys

nucs = ['A', 'C', 'G', 'T']

def targetFiles2singleCSV(dirPath, numVarPos, fileSuffix,
                          style = "verbose"):
	# Just puts a directory of nucleotide triplet binding files 
	# (NNN.txt) into one a single csv file where columns
	# are:
	# freq n1 n2 n3 a0 a1 a2 a3 a5 a6 obsCode possCode jsd entropy
	#
	# Assumes files are named <NNN.txt> and in path dir
	# Does not remove any files.
	# 
	# style = 'plain' only outputs the frequencies, nuc poisitons,
	# and amino acid positions, while 'verbose' outputs the 
	# additional information listed above.

	allFile = open(dirPath + "all.csv", 'w')
	print "Here"

	# Write the header
	if numVarPos == 6:
		if style == 'verbose':
			#allFile.write('freq,b1,b2,b3,a0,a1,a2,a3,a5,a6,' + \
	    	#          	  'obsCode,possCode,jsd,entropy\n')  # My files
			allFile.write('freq,b1,b2,b3,a0,a1,a2,a3,a5,a6,' + \
	    	          	  'possCode,obsCode,entropy\n')  #Anton files
		elif style == 'plain':
			allFile.write('freq,b1,b2,b3,a0,a1,a2,a3,a5,a6\n')
	elif numVarPos == 5:
		if style == 'verbose':
			allFile.write('freq,b1,b2,b3,a0,a2,a3,a5,a6,' + \
	    	          	'obsCode,possCode,jsd,entropy\n')
		elif style == 'plain':
			allFile.write('freq,b1,b2,b3,a0,a2,a3,a5,a6\n')

	interface = set()
	protSet = set()
	j = 0
	for n1 in nucs:
		for n2 in nucs:
			for n3 in nucs:
				triplet = n1+n2+n3
				try:
					tripletFile = open(dirPath + triplet+\
					                   fileSuffix, 'r')
				except IOError:
					print "Invalid File Name: "+triplet+\
						fileSuffix
					continue
				i = 0
				for line in tripletFile:
					sp_line = line.strip().split()
					#print sp_line
					obs = ""
					prot, freq = sp_line[0], sp_line[1]

					interface.add((triplet, prot))
					protSet.add(prot)
					
					obs += freq + ','
					for n in triplet:
						obs += n + ','
					for a in prot:
						obs += a + ','
					
					if style == 'verbose':
						#obsCode, possCode, jsd, entropy = sp_line[2], \
						#	sp_line[3], sp_line[4], sp_line[5]  #My files 
						obsCode, possCode, entropy = sp_line[3], \
							sp_line[2], sp_line[4]  #Anton files 
						#obs += ','.join([obsCode, possCode, jsd, entropy]) # My files
						obs += ','.join([possCode, obsCode, entropy]) # Anton files
						obs += '\n'
					elif style == 'plain':
						obs = obs[:-1] + '\n'
					
					allFile.write(obs)
					i += 1
				print "Num proteins in %s:  %d" \
					%(dirPath + triplet + fileSuffix, i)
				j += i
	print "Total num protein/dna observations in all files:  %d" %j
	print "Total num unique proteins in all files:  %d" %len(protSet)
	print
	allFile.close()
	return interface

def csv2txtFile(dirPath, numVarPos, style = 'verbose'):
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
			if style == 'verbose':
				#obsCode = sp_line[10]
				#possCode = sp_line[11]
				#jsd = sp_line[12]
				#entropy = sp_line[13]  # My files
				obsCode = sp_line[11]
				possCode = sp_line[10]
				entropy = sp_line[12]   # Anton files

		elif numVarPos == 5:
			prot = ''.join(sp_line[4:9])
			if style == 'verbose':
				obsCode = sp_line[9]
				possCode = sp_line[10]
				jsd = sp_line[11]
				entropy = sp_line[12]

		if style == 'verbose':
			#outStr = '\t'.join([base, prot, freq, obsCode, possCode,
		    #                jsd, entropy]) + '\n'   # My files
			outStr = '\t'.join([base, prot, freq, possCode, obsCode,
		                    	entropy]) + '\n'   # Anton files
		elif style == 'plain':
			outStr = '\t'.join([base, prot, freq]) + '\n'
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