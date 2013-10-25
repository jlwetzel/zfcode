import numpy as np
from pwm import makeLogo
import re
import os
from pwm import pwmfile2matrix, infoEntr

nucs = ['A', 'C', 'G', 'T']

def makeNucMatFile(path, label, nucMat):
	# Write a transfac style frequency matrix to file

	fout = open(path + label + '.txt', 'w')
	# Write a dummy header
	fout.write('ID idNum\nBF species\n')
	fout.write('P0\t' + '\t'.join(nucs) + '\n')
	for i in range(len(nucMat)):
		outstr = str(i+1).zfill(2)
		for j in range(len(nucMat[i,:])):
			outstr += '\t%.4f' %nucMat[i,j]
		outstr += '\tX\n'
		fout.write(outstr)
	fout.write('XX\n\\\\\n')
	fout.close()


def normalizePWM(pwm):
	# Normalizes each inner array of 
	# a 2d numpy array
	for i in range(len(pwm)):
		pwm[i,:] = pwm[i,:]/np.sum(pwm[i,:])
	return pwm

def get_GAG_pwm(pwm):
	# Find the index in the pwm where the 
	# GAG motif starts and then return a new
	# pwm that begins at this index.

	# Find the start index of the GAG motif
	G1Found = False
	A1Found = False
	GAGFound = False
	G1ind = None
	A1ind = None
	for i in range(len(pwm)):
		if GAGFound:
			break
		argmax = np.argmax(pwm[i,:])
		if G1Found and A1Found and argmax == nucs.index('G') and \
			i - 1 == A1ind:
			GAGFound = True
		elif G1Found and argmax == nucs.index('A') and \
			i - 1 == G1ind:
			A1Found = True
			A1ind = i
		elif argmax == nucs.index('G'):
			G1Found = True
			G1ind = i
	startInd = G1ind
	#print startInd

	if not GAGFound:
		return None
	if startInd == None or startInd > 3:
		return None
		
	# Copy the appropriate positions into the new PWM
	trimmedPWM = np.zeros((len(pwm) - startInd, 4), float)
	for i in range(len(trimmedPWM)):
		trimmedPWM[i,:] = pwm[i+startInd,:]

	return trimmedPWM

def get_GCG_pwm(pwm):
	# Find the index in the pwm where the 
	# GCG motif starts and then return a new
	# pwm that begins 6 positions earlier than this.

	# Find the start index of the GAG motif
	#print pwm
	G1Found = False
	C1Found = False
	GCGFound = False
	G1ind = None
	C1ind = None
	for i in range(len(pwm) - 1, 0, -1):
		if GCGFound:
			break
		argmax = np.argmax(pwm[i,:])
		if G1Found and C1Found and argmax == nucs.index('G') and \
			(i + 1) == C1ind:
			GCGFound = True
		elif G1Found and argmax == nucs.index('C') and \
			(i + 1) == G1ind:
			C1Found = True
			C1ind = i
		elif argmax == nucs.index('G'):
			G1Found = True
			G1ind = i
	
	#print G1Found
	#print C1Found
	#print GCGFound
	#print G1ind
	#print C1ind

	if not GCGFound:
		return None
	
	endInd = G1ind + 1
	startInd = endInd - 9
	if startInd == None or startInd < 0 or endInd > 11:
		return None
	
	# Copy the appropriate positions into the new PWM
	trimmedPWM = np.zeros((9, 4), float)
	for i in range(len(trimmedPWM)):
		trimmedPWM[i,:] = pwm[i+startInd,:]

	return trimmedPWM

def make3posPWM(outDir, label, pwm, targ, prot, finger):
	# Make pwms and logos for just the 3 positons
	# following the GAG motif

	pwm3 = np.zeros((3, 4), float)

	if finger == 'F2':
		for i in range(3):
			pwm3[i,:] = pwm[i+3,:]
	elif finger == 'F3':
		for i in range(3):
			pwm3[i,:] = pwm[i,:]

	makeNucMatFile(outDir + 'pwms3/', label, pwm3)
	logoIn = outDir + 'pwms3/' + label + '.txt'
	logoOut = outDir + 'logos3/' + label + '.pdf'
	makeLogo(logoIn, logoOut, alpha = 'dna',
	         format = 'pdf', 
	         colScheme = 'classic', annot = "'5,M,3'",
	         fineprint = '""')

def parseMemeFile(fpath):
	# Returns a PWM of the width of the alignment in the
	# MEME.txt file.

	print fpath
	try:
		fin = open(fpath, 'r')
	except IOError:
		print "File: %s not found!" %fpath
		return None, None

	# Find the background nuc frequencies
	for line in fin:
		sp_line = line.strip().split()
		if len(sp_line) > 1:
			if sp_line[0] == "Background":
				break				
	for i, line in enumerate(fin):
		if i > 0:
			break
		sp_line = line.strip().split()
		bgFreq = {}
		bgFreq['A'] = eval(sp_line[1])
		bgFreq['C'] = eval(sp_line[3])
		bgFreq['G'] = eval(sp_line[5])
		bgFreq['T'] = eval(sp_line[7])

	# Find the motif width
	for line in fin:
		sp_line = line.strip().split()
		if len(sp_line) > 1:
			if ' '.join([sp_line[0], sp_line[1]]) == "MOTIF 1":
				mWid = eval(sp_line[4])
				break

	# Get to the first alignment line
	countDash = 0
	for line in fin:
		line = line.strip()
		#print line
		if len(line) > 0 and line[0] == '-':
			countDash += 1
			if countDash > 5:
				break
	
	# Build the count per-position count matrix
	pwm = np.zeros((mWid, len(nucs)), float)
	for line in fin:
		if line[1] == '-':
			break
		sp_line = line.strip().split()
		if 'MNXX' not in fpath:
			count = np.log2(eval(sp_line[0].split('_')[1]))
		else:
			count = 1
		nucseq = sp_line[4]
		for i, n in enumerate(nucseq):
			pwm[i,nucs.index(n)] += count
	
	# Close the file
	fin.close()

	# Return the normalized pwm
	return normalizePWM(pwm), bgFreq

def makeallpwms(listfile, targDict, finger): 
	# Makes pwms for the barcode listfile
	# assuming we have all the information 
	# neccessary to do so
	
	# For mapping REV1, REV2, REV3 to barcode prefixes
	bcKey = {1: 'TC', 2: 'AA', 3: 'GG'}
	
	fin = open(listfile, 'r')
	fileNotFound = 0        # of times no file was found
	motif_notFound = 0      # of times no sentinel motif was found
	numPwms = 0             # total number of pwms created
	keysNotFound = []       # of times a file was found, but 
	                        # we don't have a mapping to the protein

	for i, line in enumerate(fin):
		if i > 0:
			pass#continue
		sp_line = line.strip().split()
		bc = bcKey[eval(sp_line[4][-1])] + '_' + sp_line[2]
		dset = sp_line[0].split('-')[0]
		targNum = eval(sp_line[0].split('-')[-1])
		stringency = sp_line[1]
		fpath = '../data/revExp/revExpAllData/%s/%s/' \
			%(dset, bc)

		pwm, bgFreq = parseMemeFile(fpath + 'meme.txt')
		
		# If the file was found, trim the uniformative 
		# positions from front of pwm and then write 
		# the pwm to a file in two location and create 
		# a logo.
		if pwm == None:
			fileNotFound += 1
			continue
		if finger == 'F2':
			pwm = get_GAG_pwm(pwm)
		if finger == 'F3':
			pwm = get_GCG_pwm(pwm)
		if pwm == None:
			if finger == 'F2':
				print '%smeme.txt has no valid GAG motif!' %fpath
			if finger == 'F3':
				print '%smeme.txt has no valid GCG motif!' %fpath
			motif_notFound += 1
			continue
		
		# Found a matrix with a valid GAG motif near beginning
		# Write the pwm and logo to the F2_GAG directory
		numPwms += 1
		try:
			targ = targDict[targNum][0]
			prot = targDict[targNum][1]
		except KeyError:
			keysNotFound.append(targNum)
			continue

		if finger == 'F2':
			outDir = '../data/revExp/F2_GAG/'
		elif finger == 'F3':
			outDir = '../data/revExp/F3_GCG/'

		label = '_'.join([str(targNum), targ, prot, stringency])
		makeNucMatFile(outDir + 'pwms/', label, pwm)
		logoIn = outDir + 'pwms/' + label + '.txt'
		logoOut = outDir + 'logos/' + label + '.pdf'
		makeLogo(logoIn, logoOut, alpha = 'dna',
		         colScheme = 'classic', 
		         xlab = '_'.join([targ,prot]))
		make3posPWM(outDir, label, pwm, targ, prot, finger)

	fin.close()
	print "Number of pwms created: %d" %numPwms
	print "Number of files not found: %d" %fileNotFound

	if finger == 'F2':
		print "Number where GAG motif not found: %d" %motif_notFound
	elif finger == 'F3':
		print "Number where GCG motif not found: %d" %motif_notFound

	print "Number of keys not found: %d" %len(keysNotFound)
	print "Keys not found: "
	for k in keysNotFound:
		print k
		
def getTargDict(targfname, finger):
	# Get a dictionary mapping each protNum 
	# to a (targ, protein) tuple
	targDict = {}
	fin = open(targfname, 'r')
	fin.readline()  # Skip the header
	for line in fin:
		sp_line = line.strip().split()
		targNum = eval(sp_line[0])
		if finger == 'F2':
			targ = sp_line[4].split('-')[1]
			prot = sp_line[2]
		elif finger == 'F3':
			targ = sp_line[4].split('-')[0]
			prot = sp_line[3]
		targDict[targNum] = (targ, prot)
	return targDict

def getPatternbcs(inFile, outFile, bcRegex):
	# Find the lines of inFile that match 
	# pattern regex and write only those 
	# lines to outFile

	fin = open(inFile, 'r')
	fout = open(outFile, 'w')
	for line in fin:
		sp_line = line.strip().split()
		if re.match(bcRegex, sp_line[0].split('-')[-1]) != None:
			fout.write(line)
	fout.close()
	fin.close()

def get4posF2F3RevExpInter():
	# Returns a list of fingers that are in the 
	# intersection for F2 and F3 when considering
	# only the 4 canonical positions.
	pass

def makeDir(dname):
	try:
		os.mkdir(dname)
	except OSError:
		pass

def getGoodLogoFileNames(dpath, colmin = 0.2, colAvg = 0.75):
	# Returns a list the list of files from a reverse 
	# experiment that correspond to "good" logos by some
	# measure of logo quality.  colMin is the theshold for 
	# minimum information content of each column.  colAvg
	# is the threshold for the average across the 3 columns.
	# dpath should correspond to a directory containing 
	# transfac style pwms

	handle = os.popen("ls %s" %dpath)
	keepList = []

	for fname in handle:
		fname = fname.strip()
		if '5mM' not in fname:
			continue
		nucMat = pwmfile2matrix(dpath + fname)

		keep = True
		ics = []
		for i in range(len(nucMat)):
			ic = 2 - infoEntr(nucMat[i,:])
			if ic < colmin:
				keep = False
			ics.append(ic)

		if keep and np.array(ics).mean() >= colAvg:
			keepList.append(dpath + fname)

	return keepList

def getcanHelixToFnameMap(names, canInd):

	# Place all files that map to the same 
	# canonical helix in a list keyed by the 
	# canonical helix
	hDict = {}
	for fname in names:
		prot = fname.split('/')[-1].split('_')[2]
		canHelix = ''
		for i in canInd:
			canHelix += prot[i]
		if hDict.has_key(canHelix):
			hDict[canHelix].append(fname)
		else:
			hDict[canHelix] = [fname]
	return hDict

def makeOlapPWMdir(h1, h2, sharedKeys, newDir):

	makeDir(newDir)
	pwmDir = newDir + 'pwms3/'
	logoDir = newDir + 'logos3/'
	makeDir(pwmDir)
	makeDir(logoDir)
	for k in sharedKeys:
		for fpath in h1[k]:
			fname = fpath.split('/')[-1]
			os.system("cp %s %s" \
			          %(fpath, pwmDir + k + "_F2_" + fname))
			fstem = '/'.join(fpath.split('/')[:-2]) + '/'
			pdfname = fname.split('.')[0] + '.pdf'
			os.system("cp %s %s" \
			          %(fstem + 'logos3/' + pdfname, 
			            logoDir + k + "_F2_" + pdfname))
		for fpath in h2[k]:
			fname = fpath.split('/')[-1]
			os.system("cp %s %s" \
		          	%(fpath, pwmDir + k + "_F3_" + fname))
			fstem = '/'.join(fpath.split('/')[:-2]) + '/'
			pdfname = fname.split('.')[0] + '.pdf'
			os.system("cp %s %s" \
			          %(fstem + 'logos3/' + pdfname, 
			            logoDir + k + "_F3_"+  pdfname))


def getF2F3canonHelixOlapFnames(names1, names2, canInd):
	# Get the list of finlnames for all designed 
	# helices that overlap in the canonical positions
	# (-1, 2, 3, 6)

	h1 = getcanHelixToFnameMap(names1, canInd)
	h2 = getcanHelixToFnameMap(names2, canInd)
	#print len(names1), len(names2)
	#print len(h1.keys()), len(h2.keys())
	sharedKeys = []
	for k in h1.keys():
		if k in h2.keys():
			sharedKeys.append(k)
	#print len(sharedKeys)
	#print sharedKeys

	return sharedKeys, h1, h2

def makePWMbyPosition(finger, expNum, startPos):
	# Make a 3 positions PWM from the MEME file given 
	# a starting position in the alignment

	##### NOTE ... need to fix because there can be multiple
	#####          lines for one barcode (high vs low stringency)

	bcKey = {1: 'TC', 2: 'AA', 3: 'GG'}

	# Get the right files to look stuff up
	if finger == "F2":
		bcfile = '../data/revExp/revExpBarcodes/all700Entries.txt'
		expfile = '../data/revExp/revExpBarcodes/revExper_GAG_700s.txt'
	elif finger == "F3":
		bcfile = '../data/revExp/revExpBarcodes/all100Entries.txt'
		expfile = '../data/revExp/revExpBarcodes/'

	# Look stuff up
	bcline = os.popen('grep "\\-%s" %s' %(expNum, bcfile)).readline()
	expline = os.popen('grep "^%s" %s' %(expNum, expfile)).readline()
	sp_bcline = bcline.split()
	sp_expline = expline.split()
	mn = sp_bcline[0].split('---')[0]
	tag = bcKey[eval(sp_bcline[-1][-1])] + '_' + \
		sp_bcline[-3]
	memepath = '../data/revExp/revExpAllData/%s/%s/meme.txt' \
		%(mn, tag)
	strin = sp_bcline[1]
	if finger == "F2":
		targ = sp_expline[4].split('-')[1]
		helix = sp_expline[2]
	elif finger == "F3":
		targ = sp_expline[4].split('-')[0]
		helix = sp_expline[3]

	# Parse the meme file
	pwm = parseMemeFile(memepath)[0]

	# Trim the pwm to 3 positions
	pwm3 = np.zeros((3,4), dtype = 'float')
	for i in range(len(pwm3)):
		print i
		print i + startPos
		pwm3[i,:] = pwm[(i + startPos - 1),:]

	# Write pwm to file and make logo
	label = '_'.join([expNum, targ, helix, strin])
	makeNucMatFile('./', label, pwm3)
	logoIn = label + '.txt'
	logoOut = label + '.pdf'
	makeLogo(logoIn, logoOut, alpha = 'dna',
	         format = 'pdf', 
	         colScheme = 'classic', annot = "'5,M,3'",
	         fineprint = '""')

def main():

	makePWMbyPosition('F2', '750', 6)
	"""
	# Check out how many of the logos look "good"
	F2path = '../data/revExp/F2_GAG/pwms3/'
	F3path = '../data/revExp/F3_GCG/pwms3/'
	goodF2s = getGoodLogoFileNames(F2path, 
	                               colmin = 0.2, colAvg = 0.75)
	goodF3s = getGoodLogoFileNames(F3path,
	                               colmin = 0.2, colAvg = 0.75)
	#print len(goodF2s)
	#print len(goodF3s)
	#for fname in goodF3s:
	#	print fname.split('/')[-1]

	sharedKeys, h1, h2 = getF2F3canonHelixOlapFnames(goodF2s, 
	                                                 goodF3s,
	                                     canInd = [0,2,3,6])
	newDir = "../data/revExp/F2F3sharedHelix_020_075/"
	makeOlapPWMdir(h1, h2, sharedKeys, newDir)
	"""

	"""	
	# Make the barcode files for F2 and F3 experiments
	print "Hello"
	allbcfile = '../data/revExp/revExpBarcodes/allBarcodes.txt'
	bc700sfile = '../data/revExp/revExpBarcodes/all700Entries.txt'
	bc100sfile = '../data/revExp/revExpBarcodes/all100Entries.txt'
	getPatternbcs(allbcfile, bc700sfile, r'[0-9]?7[0-9]{2}')
	getPatternbcs(allbcfile, bc100sfile, r'[0-9]?1[0-9]{2}')
	"""

	"""
	# Make the F2 pwms
	finger = 'F2'
	bcfname = '../data/revExp/revExpBarcodes/all700Entries.txt'
	targfname = '../data/revExp/revExpBarcodes/revExper_GAG_700s.txt'
	targDict = getTargDict(targfname, finger)
	#print targDict
	makeallpwms(bcfname, targDict, finger)
	"""

	"""
	# Make the F3 pwms
	finger = 'F3'
	bcfname = '../data/revExp/revExpBarcodes/all100Entries.txt'
	targfname = '../data/revExp/revExpBarcodes/revExper_GAG_100s.txt'
	targDict = getTargDict(targfname, finger)
	makeallpwms(bcfname, targDict, finger)
	"""
	


if __name__ == '__main__':
	main()