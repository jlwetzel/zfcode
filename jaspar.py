# Code for dealing with the JASPAR sql_table files and the full matrix file.

JASPAR_VERSION = '2009-Oct12'
prefix = '../data/JASPAR/' + JASPAR_VERSION + '/vertebrates/'
protTab = prefix + '/sql_tables/MATRIX_PROTEIN.txt'
annotTab = prefix + '/sql_tables/MATRIX_ANNOTATION.txt'
speciesTab = prefix + '/sql_tables/MATRIX_SPECIES.txt'
matrixTab = prefix + '/sql_tables/MATRIX.txt'
allMats = prefix + '/matrix_only.txt'

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
	return sorted(list(ids))

def JASPARIDs2proteinIDs(JASPARids):
	# Takes a sorted list of JASPAR IDs and 
	# returns a list of the corresponding protein IDs
	
	protFile = open(protTab, 'r')
	i = 0
	proteinIDs = []
	for line in protFile:
		sp_line = line.strip().split()
		if sp_line[0] == JASPARids[i]:
			proteinIDs.append(sp_line[1])
			i += 1
		if i == len(JASPARids):
			break

	protFile.close()
	return proteinIDs

def main():
	JASPARids = getIDsByAnnot(('family', 'BetaBetaAlpha-zinc finger'))
	protIDs = JASPARIDs2proteinIDs(JASPARids)
	for ID in protIDs:
		print ID

if __name__ == '__main__':
	main()