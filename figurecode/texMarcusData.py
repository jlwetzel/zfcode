import os
import re

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
				if os.path.exists(logoTypeDir[i-1] + newfname):
					texFile.write(' & ' + '\\includegraphics[scale=0.5]{%s}' \
				              	  %(logoTypeDir[i-1] + newfname))
				else:
					print "Couldn't find: " + logoTypeDir[i-1] + newfname
					texFile.write(' & ' + '- ')
		texFile.write('\\\\\n')

	texFile.write('\\end{longtable}\n')


def main():
	
	style = 'llratio'
	useExact = False
	mat = 'PAM30'
	trainFing = 'F2'
	trainStrin = 'low'
	outDir = '../../figures/predictionLogos/'
	inDirPrefix = '../../data/'

	if style == 'llratio':
		logoType = ['Exper', 'Lookup', 'NN.Top25.PAM30', 'LLRatio.Temp2']
		logoTypeDir = ['revExp/F2_GAG/logos3/', 'lookupTable/',
					   'NN_top25_' + mat + '/', 'llRatioCA/']

	if style == 'blsVnn':
		k = str(25)
		logoType = ['Exper.', 'Lookup', 'NN.Top' + k + '.' + mat, 'BLS02']
		logoTypeDir = ['revExp/F2_GAG/logos3/', 'lookupTable/', 
			       	   'NN_top' + k + '_' + mat + '/','bls02RevExp/']

	elif re.match(r'(.)*top[0-9][0-9]', style) != None and useExact:
		k = style[-2:]
		logoType = ['Exper.', 'Lookup', 'NNOnly.Top' + k + '.' + mat, \
			'NN.Top' + k + '.' + mat]
		logoTypeDir = ['revExp/F2_GAG/logos3/', 'lookupTable/', 
			       	   'NNonly_top' + k + '_' + mat + '/','NN_top' + k + '_' + mat + '/']

	elif re.match(r'(.)*top[0-9][0-9]', style) != None:
		k = style[-2:]
		logoType = ['Exper.', 'Lookup', 'NNOnly.Top' + k, \
			'NNOnly.Top' + k + '.' + mat]
		logoTypeDir = ['revExp/F2_GAG/logos3/', 'lookupTable/', 
			       	   'NNonly_top' + k + '/','NN_top' + k + '_' + mat + '/']

	#filters = ['cut3bc_0_5', 'cut10bc_0_5', 'cut10bc_0', 'cut10bc_025', 'cut3bc_025']
	filters = ['cut10bc_0_5']
	for f in filters:
		if useExact:
			texFile = open(outDir + '_'.join([trainFing, trainStrin, style, f, 'useExact']) \
		         	       + '.tex', 'w')
		else:
			texFile = open(outDir + '_'.join([trainFing, trainStrin, style, f]) \
		    	           + '.tex', 'w')
		makeTexPreamble(texFile)
		logoTypeDir2 = logoTypeDir[:]
		for i, l in enumerate(logoTypeDir):
			logoTypeDir2[i] = inDirPrefix + logoTypeDir2[i]
			if i != 0:
				if style == 'blsVnn' and i == 3:
					logoTypeDir2[i] += '/'.join(['F2', 'low', \
					                            'predictions', 'logos']) + '/'
				else:
					logoTypeDir2[i] += '/'.join([trainFing, trainStrin, f,
				    	                        'predictions', 'logos']) + '/'
		
		#print logoTypeDir2
		makeTexTable(texFile, logoType, logoTypeDir2)
		makeTexEnding(texFile)

		texFile.close()

if __name__ == '__main__':
	main()