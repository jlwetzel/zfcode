# Script for manipulating parsing through the B1H database
# starting at the level of variable protein regions in the 
# form of codons.

import os

# Global Constants
codon2amino = {'TTT':  'F', 'TTC':  'F', 'TTA':  'L', 'TTG':  'L',
	           'TCT':  'S', 'TCC':  'S', 'TCA':  'S', 'TCG':  'S', 
	           'TAT':  'Y', 'TAC':  'Y', 'TAA':  '*', 'TAG':  '*', 
	           'TGT':  'C', 'TGC':  'C', 'TGA':  '*', 'TGG':  'W', 
	           'CTT':  'L', 'CTC':  'L', 'CTA':  'L', 'CTG':  'L', 
	           'CCT':  'P', 'CCC':  'P', 'CCA':  'P', 'CCG':  'P',
	           'CAT':  'H', 'CAC':  'H', 'CAA':  'Q', 'CAG':  'Q', 
	           'CGT':  'R', 'CGC':  'R', 'CGA':  'R', 'CGG':  'R', 
	           'ATT':  'I', 'ATC':  'I', 'ATA':  'I', 'ATG':  'M', 
	           'ACT':  'T', 'ACC':  'T', 'ACA':  'T', 'ACG':  'T', 
	           'AAT':  'N', 'AAC':  'N', 'AAA':  'K', 'AAG':  'K', 
	           'AGT':  'S', 'AGC':  'S', 'AGA':  'R', 'AGG':  'R', 
	           'GTT':  'V', 'GTC':  'V', 'GTA':  'V', 'GTG':  'V',
	           'GCT':  'A', 'GCC':  'A', 'GCA':  'A', 'GCG':  'A', 
	           'GAT':  'D', 'GAC':  'D', 'GAA':  'E', 'GAG':  'E',
	           'GGT':  'G', 'GGC':  'G', 'GGA':  'G', 'GGG':  'G'}

codonBias = {'A':  {'GCT':.19, 'GCC':.25, 'GCA':.22 , 'GCG':.34},
             'C':  {'TGT':.43, 'TGC':.57},
             'D':  {'GAT':.59, 'GAC':.41},
             'E':  {'GAA':.7, 'GAG':.3},
             'F':  {'TTT':.51, 'TTC':.49},
             'G':  {'GGT':.38, 'GGC':.4, 'GGA':.09, 'GGG':.13},
             'I':  {'ATT':.47, 'ATC':.46, 'ATA':.07},
             'H':  {'CAT':.52, 'CAC':.48},
             'K':  {'AAA':.76, 'AAG':.24},
             'L':  {'TTA':.11, 'TTG':.11, 'CTT':.1, 'CTC':.1, 
                    'CTA':.03, 'CTG':.55},
             'M':  {'ATG':1.0},
             'N':  {'AAT':.39, 'AAC':.61},
             'P':  {'CCT':.16, 'CCC':.1, 'CCA':.2, 'CCG':.55},
             'Q':  {'CAA':.31, 'CAG':.69},
             'R':  {'CGT':.42, 'CGC':.37, 'CGA':.05, 'CGG':.08, 
                    'AGA':.04, 'AGG':.03},
             'S':  {'TCT':.19, 'TCC':.17, 'TCA':.12, 'TCG':.13, 
                    'AGT':.13, 'AGC':.27},
             'T':  {'ACT':.18, 'ACC':.37, 'ACA':.26, 'ACG':.2},
             'V':  {'GTT':.29, 'GTC':.2, 'GTA':.17, 'GTG':.34},
             'W':  {'TGG':1.0},
             'Y':  {'TAT':.53, 'TAC':.47},
             '*':  {'TAA':.62, 'TAG':.09, 'TGA':.3}}

class TargetObj(object):
	def __init__(self, targ, valList):
		# targ is the target 3-mer
		# valList is a list of lists of the form:
		# [[seqrun, barcode, finger, stringency], ... ]
		self.targ = targ
		self.seqruns, self.bcodes, self.fings, self.strins = \
			[], [], [], []
		for l in valList:
			self.seqruns.append(l[0]) 
			self.bcodes.append(l[1])
			self.fings.append(l[2])
			self.strins.append(l[3])
	
	def getValList(self):
		# Returns a list of lists of the same form 
		# that was passed into constructor
		return zip(self.seqruns, self.bcodes, 
		           self.fings, self.strins)

	def addValsToList(self, valList):
		# For adding a new experiment to the target object
		# valList is a lists of the form:
		# [seqrun, barcode, finger, stringency]
		self.seqruns.append(l[0]) 
		self.bcodes.append(l[1])
		self.fings.append(l[2])
		self.strins.append(l[3])

	def getByFingAndStrin(self, fing, strin):
		# Returns a list of seqruns that are associated
		# with the finger and stringency for this target
		seqruns = []
		for i in range(len(self.seqruns)):
			if self.fings[i] == fing and \
				self.strins[i] == strin:
				seqruns.append(self.seqruns[i])
		return seqruns
			
def parseMasterFile(path):
	# Returns a dictionary of 3merTargetObjs indexed by 
	# their 3mer target strings
	mfile = open(path, 'r')
	
	targs = {}
	# Create a dictionary: target -> list of 4-tuples
	for line in mfile:
		if line == '' or line == '\n':
			continue
		elif line[0] == '>':
			sp_line = line.strip().split()
			finger = sp_line[0][1:]
			stringency = sp_line[1]
		else:
			sp_line = line.strip().split('\t')
			if len(sp_line) > 1:
				target = sp_line[0]
				for i in range(1, len(sp_line)):
					seqrun, barcode = sp_line[i].strip().split()
					if targs.has_key(target):
						targs[target].append([seqrun, barcode,
						                      finger, stringency])
					else:
						targs[target] = [[seqrun, barcode,
						                  finger, stringency]]

	# Create the dict of TargObjs
	targObjs = {}
	for t in targs.keys():
		targObjs[t] = TargetObj(t, targs[t])

	# Close file and return dict of TargObjs
	mfile.close()
	return targObjs

def organizeByTarget(targObjs, newpath, oldpath, mode = 'c'):
	# Creates a new directory structure for the 
	# NNN_all_nucleotide_sequences.txt files with 
	# a hierarchy of:
	# finger -> stringency -> NNN_MNXX_all_nuc_seq.txt
	#
	# mode is 'c' by default to copy the files. 
	# set mode to 'm' if you want to move the files
	# instead of copying.
	#
	# The root of the new dir structure is given 
	# by newpath.  The root of the old dir structure
	# is given by oldpath.
	
	try:
		os.mkdir(newpath)
	except OSError:
		pass
	
	fings = ['F1', 'F2', 'F3']
	strins = ['low', 'high']

	# Make the directory hierarchy	
	for fing in fings:
		try:
			os.mkdir(newpath + '/' + fing + '/')
		except OSError:
			pass

		for strin in strins:
			try:
				os.mkdir(newpath + '/' + fing + '/' + \
				         strin + '/')
			except OSError:
				pass

			# Copy the files for each target
			for targ in sorted(targObjs.keys()):

				seqruns = targObjs[targ].getByFingAndStrin(fing, strin)
				for seqrun in sorted(seqruns):
					oldfile = oldpath + '/' + seqrun + '/' + \
				         	  fing + '/' + targ \
				         	  + '_all_nucleotide_sequences.txt'
					newfile = newpath + '/' + fing + '/' + \
				         	  strin + '/'+ targ + '_' + seqrun \
				         	  + '_all_nuc_seq.txt'
					if mode == 'c':
						print "Copying %s to %s." \
							%(oldfile, newfile)
						os.system('cp ' + oldfile + ' ' + newfile)
					elif mode == 'm':
						print "Moving %s to %s." \
							%(oldfile, newfile)
						os.system('mv ' + oldfile + ' ' + newfile)

				# Correct the permissions for the files
				newdir = newpath + '/' + fing + '/' + strin + '/'
				os.system('chmod 644 ' + newdir + '*')

def main():
	mFilePath = '../data/b1hdata/database/all_selection_data_final.txt'
	
	# Get the dict of TargObjs
	targs = parseMasterFile(mFilePath)
	
	# Test the class
	#for k in sorted(targs.keys()):
	#	print "Targ:  %s" %k
	#	print targs[k].getByFingAndStrin('F3', 'low')

	# Make the new directory structure
	oldpath = "../data/b1hData/database"
	newpath = "../data/b1hData/newDatabase/6varpos"
	organizeByTarget(targs, newpath, oldpath, mode = 'c')

if __name__ == '__main__':
	main()