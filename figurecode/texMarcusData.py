import os

def makeTexPreamble(texFile):
	preamble = '\n'.join(['\\documentclass[11pt]{article}',
						  '\\usepackage[latin1]{inputenc}',
						  '\\usepackage{amsmath}',
						  '\\usepackage{amsfonts}',
						  '\\usepackage{amssymb}',
						  '\\usepackage{verbatim}',
						  '\\usepackage{geometry}',
						  '\\usepackage{booktabs}',
						  '\\usepackage{longtable}',
						  '\\usepackage{color}',
						  '\\usepackage{graphicx}']) + '\n\n'
	texFile.write(preamble)
	texFile.write('\\begin{document}\n\n')

def makeTexEnding(texFile):
	texFile.write('\\end{document}\n')

def makeTexTable(texFile, logoType, logoTypeDir):
	# Takes all of the logos from a 
	# directory and puts them all into one 
	# large latex table

	colStr = ''
	for i in range(len(logoType) + 1):
		colStr += 'c'

	texFile.write('\\begin{longtable}{%s}\n' %(colStr))

	ncol = len(logoType)+ 1
	flist = os.popen('ls ' + logoTypeDir[0])
	flist = [i.strip() for i in flist if '5mM' in i]

	# Make the column header
	texFile.write('\\hline\n')
	for i in range(ncol):
		if i == 0:
			pass
		else:
			texFile.write(' & %s' %logoType[i - 1])
	texFile.write('\\\\\n')
	texFile.write('\\hline\n')

	# Make the columns
	for fname in flist:
		logoNum = fname.split('_')[0]
		targ = fname.split('_')[1]
		prot = fname.split('_')[2]
		strin = fname.split('_')[3].split('.')[0]

		for i in range(ncol):
			if i == 0:
				texFile.write('.'.join([logoNum, strin]))
			elif i == 1:
				texFile.write(' & ' + '\\includegraphics[scale=0.5]{%s}' \
				              	  %(logoTypeDir[i-1] + fname))
			else:
				newfname = '_'.join([logoNum, targ, prot, 'low']) + '.pdf'
				print logoTypeDir[i-1] + newfname
				if os.path.exists(logoTypeDir[i-1] + newfname):
					texFile.write(' & ' + '\\includegraphics[scale=0.5]{%s}' \
				              	  %(logoTypeDir[i-1] + newfname))
				else:
					texFile.write(' & ' + '- ')
		texFile.write('\\\\\n')

	texFile.write('\\end{longtable}\n')


def main():
	outDir = '../../figures/predictionLogos/'
	inDirPrefix = '../../data/'
	logoType = ['exp', 'lookup', 'lookNNTripP30', 'N.lookNNTrip', \
				'N.lookNNTripP30']
	logoTypeDir = ['revExp/F2_GAG/logos3/', 'lookupTable/', 
			       'lookupTableNNonly_pam30_decomp2/',
			       'lookNNonly_Triples/', 'lookNNonly_PAM30_Triples/']
	fing = 'F2'
	strin = 'low'
	filters = ['cut3bc_0_5']#['cut10bc_0_5, cut10bc_0']
	for f in filters:
		texFile = open(outDir + '_'.join([fing, strin, f]) + '.tex', 'w')
		
		makeTexPreamble(texFile)
		logoTypeDir2 = logoTypeDir[:]
		for i, l in enumerate(logoTypeDir):
			logoTypeDir2[i] = inDirPrefix + logoTypeDir2[i]
			if i != 0:
				logoTypeDir2[i] += '/'.join([fing, strin, f,
				                            'predictions', 'logos']) + '/'
		
		makeTexTable(texFile, logoType, logoTypeDir2)
		makeTexEnding(texFile)

		texFile.close()


if __name__ == '__main__':
	main()