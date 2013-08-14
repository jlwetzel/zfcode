import os
import sys
import re

def makeTexTable(logoDir, f, s, t):
	# Takes all of the logos from a 
	# directory and puts them all into one 
	# large latex table

	try:
		texFile = open(logoDir + 'allLogos.tex', 'w')
	except IOError:
		os.mkdir(logoDir)
		texFile = open(logoDir + 'allLogos.tex', 'w')

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
	texFile.write('%s %s %s\n' %(f,s,t.split('_')[-1]))
	texFile.write('\\begin{longtable}{ccc}\n')

	ncol = 3
	handle = os.popen('ls ' + logoDir)
	i = 0
	for fname in handle:
		fname = fname.strip()
		if re.match(r'[ACGT]{3}(.)*.pdf', fname) == None:
			continue
		texFile.write('\\includegraphics[scale=0.9]{%s}' \
		              %fname)
		if i % ncol == ncol - 1:
			texFile.write('\\\\\n')
		else:
			texFile.write('&\n')
		i += 1

	texFile.write('\\end{longtable}\n')
	texFile.write('\\end{document}\n')

	texFile.close()
	
def main():
	prefix = '../../data/b1hData/newDatabase/6varpos/'
	fings = ['F2']
	strins = ['low', 'high']
	thresh = ['protein_seq_cut3bc_0_5']
	for f in fings:
		for s in strins:
			for t in thresh:
				logoDir = prefix + '/'.join([f,s,t]) + '/logos/'
				print logoDir
				makeTexTable(logoDir, f, s, t)

if __name__ == '__main__':
	main()
