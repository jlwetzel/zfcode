# Code for reteiving and maniupulating the JASPAR sql_table files 
# and the JASPAR PWM file.

import os

JASPAR_BUILD = '2009-Oct12-NonRedundant'

prefix = '../data/JASPAR/' + JASPAR_BUILD
protTab = prefix + '/sql_tables/MATRIX_PROTEIN.txt'
annotTab = prefix + '/sql_tables/MATRIX_ANNOTATION.txt'
speciesTab = prefix + '/sql_tables/MATRIX_SPECIES.txt'
matrixTab = prefix + '/sql_tables/MATRIX.txt'
PWMfile = prefix + '/matrix_only.txt'

def getNewBuild():
	# Get the latest build of the complete JASPAR CORE set.
	# First set up directory structure in ../data/JASPAR/
	JASPAR_HTML_PREFIX = "http://jaspar.genereg.net//" + \
		"html/DOWNLOAD/jaspar_CORE/non_redundant/all_species/"
	sqlTables = ["MATRIX.txt", "MATRIX_ANNOTATION.txt", "MATRIX_DATA.txt",
			 	"MATRIX_PROTEIN.txt", "MATRIX_SPECIES.txt"]

	os.mkdir("../data/JASPAR/" + JASPAR_BUILD)
	os.mkdir("../data/JASPAR/" + JASPAR_BUILD + "/sql_tables")
	for tab in sqlTables:
		os.system("wget -P " + prefix + "/sql_tables/ " + 
		          JASPAR_HTML_PREFIX + "/sql_tables/" + tab) 
	os.system("wget -P " + prefix + " " + JASPAR_HTML_PREFIX 
	          + "matrix_only/matrix_only.txt") 


def getIDsByAnnot(annot, currentList = None):
	# Returns a list of JASPAR unique IDs that are are
	# labelled by the annots.  annots is tuple (key, value)
	if currentList == None:
		ids = set()
	else:
		ids = set(currentList)

	annotFile = open(annotTab, 'r')
	for line in annotFile:
		sp_line = line.strip().split('\t')

		if len(sp_line) < 3:
			continue
		key = sp_line[1]
		val = sp_line[2]
		if key == annot[0] and val == annot[1]:
			ids.add(sp_line[0])

	annotFile.close()
	ids = list(ids)
	ids = [int(i) for i in ids]
	return sorted(list(ids))

def JASPARIDs2proteinIDs(JASPARids):
	# Takes a sorted list of JASPAR IDs and 
	# returns a list of the corresponding protein IDs
	
	protFile = open(protTab, 'r')
	i = 0
	proteinIDs = []
	for line in protFile:
		sp_line = line.strip().split()
		if int(sp_line[0]) == JASPARids[i]:
			proteinIDs.append(sp_line[1])
			i += 1
		if i == len(JASPARids):
			break

	protFile.close()
	return proteinIDs

def getAnnotsByJASPARid(JASPARids, label):
	# Finds the annotation associated with the JasparID
	# and label for each ID in the ***SORTED*** 
	# list of sorted JASPARids	
	annotFile = open(annotTab, 'r')
	i = 0
	vals = []
	
	for line in annotFile:
		if len(line) != 0:
			sp_line = line.strip().split('\t')
		if int(sp_line[0]) > JASPARids[i]:
			print "No label: %s for JASPAR id %d" %(label, JASPARids[i])
			i += 1
			if i == len(JASPARids):
				break
		if int(sp_line[0]) == JASPARids[i] and sp_line[1] == label:
			vals.append(sp_line[2])
			i += 1
		if i == len(JASPARids):
			break

	annotFile.close()
	return vals

def main():
	#getNewBuild()
	JASPARids = getIDsByAnnot(('family', 'BetaBetaAlpha-zinc finger'))
	print JASPARids
	x = getAnnotsByJASPARid(JASPARids, "family")
	#protIDs = JASPARIDs2proteinIDs(JASPARids)
	#print(len(protIDs))
	for t in x:
		print t

if __name__ == '__main__':
	main()