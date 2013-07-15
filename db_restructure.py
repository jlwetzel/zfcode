# Script for processing the restructuring the database
# from the old hierarchy to the new one.
# old:  db -> seqrun -> finger -> barcode labelled sequence files.
# new:  db -> varpos -> finger -> stringency -> target labelled 
#                                               seq files.
# Also contains a function for converting a master file of 
# barcodes info to a nicely structured table.

import os

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
		
		self.bcodeDict = {}
		for i in range(len(self.seqruns)):
			if self.bcodeDict.has_key((self.seqruns[i], self.fings[i],
			                		   self.strins[i])):
				self.bcodeDict[(self.seqruns[i], self.fings[i],
			                    self.strins[i])].append(self.bcodes[i])
			else:
				self.bcodeDict[(self.seqruns[i], self.fings[i],
			                    self.strins[i])] = [self.bcodes[i]]
	
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

	# For each finger	
	for fing in fings:
		try:
			os.mkdir(newpath + '/' + fing + '/')
		except OSError:
			pass

		# For each stringency
		for strin in strins:
			try:
				os.mkdir(newpath + '/' + fing + '/' + \
				         strin + '/')
			except OSError:
				pass
			try:
				os.mkdir(newpath + '/' + fing + '/' + \
				         strin + '/all_nuc_seq/')
			except:
				pass

			# For each 3-mer target
			for targ in sorted(targObjs.keys()):

				# For each sequencing run for this 3-mer target
				seqruns = targObjs[targ].getByFingAndStrin(fing, strin)
				for seqrun in sorted(seqruns):

					# Get the correct barcodes for this particular
					# combination of target, stringency, finger,
					# and sequencing run.
					bcodes = targObjs[targ].bcodeDict[(seqrun, fing,
					                                   strin)]

					for bcode in bcodes:
						# Construct path (labelled by barcode)
						oldfile = oldpath + '/' + seqrun + '/' + \
				         	  	fing + '/' + bcode \
				         	  	+ '_all_nucleotide_sequences.txt'

				    	# Construct the new path (labelled by 3-mer target)
						newfile = newpath + '/' + fing + '/' + \
				         	  	strin + '/all_nuc_seq/'+ targ + '_' \
				         	  	+ seqrun + '_bc' + bcode + \
				         	  	'_all_nuc_seq.txt'

				    	# Copy or move the old file to new location
						if mode == 'c':
							print "Copying %s to %s." \
								%(oldfile, newfile)
							os.system('cp ' + oldfile + ' ' + newfile)
						elif mode == 'm':
							print "Moving %s to %s." \
								%(oldfile, newfile)
							os.system('mv ' + oldfile + ' ' + newfile)

				# Correct the permissions for the files
				newdir = newpath + '/' + fing + '/' + strin + \
					'/all_nuc_seq/'
				os.system('chmod 644 ' + newdir + '*')

def main():
	
	# Test the class
	#for k in sorted(targs.keys()):
	#	print "Targ:  %s" %k
	#	print targs[k].getByFingAndStrin('F3', 'low')

	# Restructure the db for the 6 variable position data
	# Get the dict of TargObjs
	mFilePath = '../data/b1hdata/database/' +\
		'all_selection_data_final.txt'
	targs = parseMasterFile(mFilePath)
	oldpath = "../data/b1hData/database"
	newpath = "../data/b1hData/newDatabase/6varpos"
	organizeByTarget(targs, newpath, oldpath, mode = 'c')

	"""
	# Restructure the db for the 5 variable position data
	# Get the dict of TargObjs
	mFilePath = '../data/b1hdata/database/' +\
		'all_selection_data_final_5pos.txt'
	targs = parseMasterFile(mFilePath)
	oldpath = "../data/b1hData/database"
	newpath = "../data/b1hData/newDatabase/5varpos"
	organizeByTarget(targs, newpath, oldpath, mode = 'c')
    """

if __name__ == '__main__':
	main()