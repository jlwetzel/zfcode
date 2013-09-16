# A module for converting cleaned protein data files into the 
# format for using SVM Light for regression.

nucs = ['A', 'C', 'G', 'T']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
	      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def getContactPairMap(contacts, baseInputMap, aminoInputMap):
	"""
	Returns a mapping of string index pairs for base-position
	and contact position using the format:

	((baseIndex, base), (helixIndex, amino)): uniqueNum

	- contacts is the name of the structural model (see makeSVMLightInput)
	- baseInputMap maps base positions to indices in the input file
	- aminoInputMap maps helix positions to indicies from the input file
	"""

	cPairMap {}
	if contacts == "ca":
		pass



def makeSVMLightInput(inPath, outPath, baseInputMap, contacts = "ca"):
	"""
	Converts a the protein input file located at inPath to the
	format recognized by SVM Lite and places it in the file 
	named outPath.  

	- Contacts is the name given to the structural model to be 
	  considered when deriving features.
	  "c" is only the canonical contacts
	  "ca" is canonical contacts + Anton's contact
	  "allPairs" is all pairwise contacts available
	"""

	contactMap = getContactPairMap(contacts, baseInputMap, aminoInputMap)

def main():
	makeSVMLightInput("../data/b1hData/newDatabase/6varpos/F2/" + \
	                  "low/protein_seq_cut10bc_025/")
