# A module for converting cleaned protein data files into the 
# format for using SVM Light for regression.


nucs = ['A', 'C', 'G', 'T']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
	      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Maps base or helix positions (0 corresponds to -1 for helix)
# to indicies of strings from the cleaned input files
baseposMap = {1:0, 2:1, 3:2}
aminoposMap = {0:0, 1:1, 2:2, 3:3, 5:4, 6:5}

def getContactPairMap(contacts):
	"""
	Returns a mapping of string index pairs for base-position
	and contact position using the format:

	(basepos, helixpos, base, amino): varName

	- contacts is the name of the structural model (see makeSVMLightInput)
	- baseInputMap maps base positions to indices in the input file
	- aminoInputMap maps helix positions to indicies from the input file
	"""

	cpairMap = {}

	if contacts == 'c':
		pairs = [(1,6), (2,3), (3,0)]
	elif contacts == 'ca':
		pairs = [(1,6), (2,3), (3,2), (3,0)]
	elif contacts == 'all':
		pairs = []
		for bpos in sorted(baseposMap.values()):
			for apos in sorted(aminoposMap.values()):
				pairs.append((bpos, apos))

	for (bpos, apos) in pairs:
		for b in nucs:
			for a in aminos:
				varName = 'b%da%d%s%s' %(bpos, apos, b, a)
				cpairMap[(bpos, apos, b, a)] = varName

	return cpairMap

def lineToBinary(line, cpairMap):
	"""
	Turn the line into an obsevation for the matrix.
	"""

	sp_line = line.strip().split()
	dna = sp_line[0]
	prot = sp_line[1]
	freq = sp_line[2]

	obs = freq
	for (bpos, apos, b, a) in sorted(cpairMap.keys()):
		if dna[baseposMap[bpos]] == b and prot[aminoposMap[apos]] == a:
			obs += ',1'
		else:
			obs += ',0'
	obs += '\n'

	return obs

def binaryToContactList(binLine, cpairMap):
	"""
	Convert from the binary line representation back to 
	the list of contatcs observed for each position in 
	the contact pair map.
	"""

	sp_line = binLine.strip().split(',')
	freq = eval(sp_line[0])
	sp_line = sp_line[1:]
	
	contactList = []
	for i, k in enumerate(sorted(cpairMap.keys())):
		if sp_line[i] == "1":
			contactList.append(cpairMap[k])

	return freq, tuple(contactList)

def makeMatrix(inPath, outPath, contacts = "ca",
               combine = "none"):
	"""
	Converts a the protein input file located at inPath to the
	format recognized for a corresponding R dataframe and places 
	it in the file named outPath.  

	- Contacts is the name given to the structural model to be 
	  considered when deriving features.
	  "c" is only the canonical contacts
	  "ca" is canonical contacts + Anton's contact
	  "allPairs" is all pairwise contacts available
	- Combine refers to how observations with the same feature 
	  vector should be combined.  Default is None, where
	  each observation is left unchanged.  
	  Also can use: "add" in which case the frequency for all 
	  lines from with identical feature vectors are added together.
	"""

	cpairMap = getContactPairMap(contacts)
	fin = open(inPath, 'r')
	fout = open(outPath, 'w')
	
	# Make the header
	headStr = 'freq'
	for k in sorted(cpairMap.keys()):
		headStr += ',' + cpairMap[k]
	headStr += '\n'
	fout.write(headStr)

	obsDict = {}
	for i, line in enumerate(fin):

		if combine == None:
			# Assemble observation in terms of binary variables
			obs = lineToBinary(line, cpairMap)
			#print i, binaryToContactList(obs, cpairMap)
			fout.write(obs)

		# Combine frequencies of identical feature vectors together
		# into a dictionary of 
		elif combine == "add":
			obs = lineToBinary(line, cpairMap)
			featureStr = ','.join(obs.strip().split(',')[1:])
			freq = eval(obs.strip().split(',')[0])
			if obsDict.has_key(featureStr):
				obsDict[featureStr] += freq
			else:
				obsDict[featureStr] = freq


	if combine == "add":
		for k in obsDict.keys():
			obs = str(obsDict[k]) + ',' + k + '\n'
			fout.write(obs)

	fout.close()
	fin.close()

def main():
	inPref = "../data/b1hData/newDatabase/6varpos/F2/" + \
	           "low/protein_seq_cut10bc_025/"
	makeMatrix(inPref + 'all.txt', inPref + 'all_matrix_ca_add.csv', 
	           contacts = 'ca', combine = "add")

if __name__ == '__main__':
	main()
