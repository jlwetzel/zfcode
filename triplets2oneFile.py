# Just puts a directory of nucleotide triplet binding files 
# (NNN.txt) into one a single csv file where columns
# are:
# count n1 n2 n3 a0 a1 a2 a3 a5 a6
#
# Assumes files are named <NNN.txt> and in path dir

import os
import sys

nucs = ['A', 'C', 'G', 'T']
path = sys.argv[1]

def main():
	os.chdir(path)
	allFile = open("all.csv", 'w')
	allFile.write("count,n1,n2,n3,a0,a1,a2,a3,a5,a6\n")
	for n1 in nucs:
		for n2 in nucs:
			for n3 in nucs:
				triplet = n1+n2+n3
				try:
					tripletFile = open(triplet+'.txt', 'r')
				except IOError:
					print "Invalid File Name: "+triplet+'.txt'
					pass
				for line in tripletFile:
					obs = ""
					(prot, count) = line.strip().split()
					obs += count + ","
					for n in triplet:
						obs += n + ","
					for a in prot[:-1]:
						obs += a + ","
					obs += a[-1] + "\n"
					allFile.write(obs)
	allFile.close()

main()