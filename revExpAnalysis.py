# A module for predicting binding pwms for the
# proteins we have reverse experiments for.

from lookupTable2 import *

def lookupMarcusPWMs(inDir, outputDir, freqDict,
                     finger, strin, filt, pred, useNN = True,
                     skipExact = False, decompose = None,
                     topk = None):
	
	# Make predcitions for each of the proteins that 
	# Marcus made experimental PWMs for

	# Create the prediction directory structure
	predictionDir = outputDir+'predictions/'
	pwmdir = predictionDir + 'pwms/'
	logodir = predictionDir + 'logos/'
	makeDir(predictionDir)
	makeDir(pwmdir)
	makeDir(logodir)
	
	# Get the directory containing the 3 position
	# experimental PWMs for this finger and set up
	# the output file for writing results
	if finger == 'F2':
		expDir = '../data/revExp/F2_GAG/pwms3/'
	elif finger == 'F3':
		expDir = '../data/revExp/F3_GCG/pwms3/'
	#print predictionDir
	fout = open(predictionDir + 'compare.txt', 'w')
	# Write header to results file
	fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
	           %('num', 'targ','prot','canonprot','c1pcc','c2pcc',
	             'c3pcc', 'c1pcc.ic', 'c2pcc.ic', 'c3pcc.ic', 'p.c1cons',
	             'p.c2cons', 'p.c3cons', 'e.c1cons', 'e.c2cons', 'e.c3cons',
	             'pred.filt', 'finger', 'strin'))
	
	
	# B1H forward experiments.  Need to update if start 
	# using the 5 position ones
	npos = 6
	canonical = True
	ind = getPosIndex(npos, canonical)
	
	# Try to do a lookup table prediction for each of the 
	# proteins for which we have an experimental 3pos pwm
	for fname in os.popen('ls ' + expDir):

		# Skip if this is the wrong stringency.
		fname = fname.strip()
		
		#if fname != "2702_AAC_EATSLRN_5mM.txt":
		#	continue


		if strin == 'low' and re.match(r'(.)*_5mM.txt', fname) == None:
			continue
		if strin == 'high' and re.match(r'(.)*_20mM.txt', fname) == None:
			continue

		# Get the info about this experiment
		sp_fname = fname.split('_')
		protNum = sp_fname[0]
		goal = sp_fname[1]
		prot = sp_fname[2].split('.')[0]
		canonProt = prot[0] + prot[2] + prot[3] + prot[6]
		label = '_'.join([str(protNum), goal, prot, strin])
	
		targList, neighborPerBase = get3merList(freqDict, canonProt, canonical,
		    	               					useNN, skipExact, decompose,
		    	               					topk)
		
		# This is target/frequncy list
		if isinstance(targList, list):
			# If the target list is empty don't output anything
			if targList == []:
				continue
			
			# Convert the 3mer target list to a frequency matrix,
			# write to file, and create a logo
			nucMat = targListToFreqMat(targList)
		
		# This is a proper numpy array already
		else:
			nucMat = targList
			#print nucMat

		makeNucMatFile(pwmdir, label, nucMat)
		logoIn = pwmdir + label + '.txt'
		logoOut = logodir + label + '.pdf'
		makeLogo(logoIn, logoOut, alpha = 'dna', 
		         colScheme = 'classic',
		         annot = "'5,M,3'",
		         xlab = '_'.join([goal,prot]))
		
		# Compare this pwm to the reverse experiment
		expMat = pwmfile2matrix(expDir + fname)
		colPcc, colPcc_ic = comparePCC(nucMat, expMat)
		predCons = getConsensus(nucMat)
		#print predCons
		expCons = getConsensus(expMat)

		# Write the comparison results to file
		fout.write("%s\t%s\t%s\t%s\t" %(protNum, goal, prot, canonProt))
		fout.write("%.4f\t"*6 %(colPcc[0], colPcc[1], colPcc[2], \
		           				colPcc_ic[0], colPcc_ic[1], colPcc_ic[2]))
		fout.write("%s\t"*6 %(predCons[0], predCons[1], predCons[2],
		           			  expCons[0], expCons[1], expCons[2]))
		#print pred+'.'+filt
		fout.write("%s\t%s\t%s\n" %(pred+'.'+filt, finger, strin))

	fout.close()

def getLabels(style, decomp, weight_mat, useExact):
	# Returrns a label and a directory prefix for the 
	# given style, decomposition, and weight matrix
	if style == 'nnonly' and useExact:

		if decomp == 'triples' and weight_mat == 'PAM30':
			label = 'NNOnly.trip.PAM30'
			inDirPref = '../data/NNonly_Triples_PAM30/'
		elif decomp == 'triples':
			label = 'NNOnly.trip'
			inDirPref = '../data/NNonly_Triples/'
		elif decomp == 'doubles' and weight_mat == 'PAM30':
			label = 'NNOnly.doub.PAM30'
			inDirPref = '../data/NNonly_Doubles_PAM30/'
		elif decomp == 'doubles':
			label = 'NNOnly.doub'
			inDirPref = '../data/NNonly_Doubles/'
		elif decomp == 'singles' and weight_mat == 'PAM30':
			label = 'NNOnly.sing.PAM30'
			inDirPref = '../data/NNonly_Singles_PAM30/'
		elif decomp == 'singles':
			label = 'NNOnly.sing'
			inDirPref = '../data/NNonly_Singles/'

	elif re.match(r'top[0-9][0-9]', style) != None and useExact:
		if weight_mat == 'PAM30':
			label = 'NN.' + style + '.PAM30'
			inDirPref = '../data/NN_' + style + '_PAM30/'
		else:
			label = 'NN.' + style 
			inDirPref = '../data/NN_' + style + '/'

	elif re.match(r'top[0-9][0-9]', style) != None:
		if weight_mat == 'PAM30':
			label = 'NNonly.' + style + '.PAM30'
			inDirPref = '../data/NNonly_' + style + '_PAM30/'
		else:
			label = 'NNonly.' + style 
			inDirPref = '../data/NNonly_' + style + '/'
		#print weight_mat
		#print label

	elif style == 'lookonly':
		label = 'look'
		inDirPref = '../data/lookupTable/'
				
	return label, inDirPref

def getDecompDict(style, decomp):
	# Return the decomposition dictionary that corresponds 
	# to the string decomp

	triples = {1: [1,2,3], 2: [1,2,3], 3: [0,1,2]}
	doubles = {1: [2,3], 2: [2,3], 3: [0,1]}
	singles = {1: [3], 2: [2], 3: [0]}

	if re.match(r'top[0-9][0-9]', style) != None:
		return singles

	elif decomp == "triples":
		return triples
	elif decomp == "doubles":
		return doubles
	elif decomp == "singles":
		return singles
	
	else:
		return None

def runBLS02Analysis(outputDir, strin):
	# Make predcitions for each of the proteins that 
	# Marcus made experimental PWMs for

	# Create the prediction directory structure
	predictionDir = outputDir+'predictions/'
	pwmdir = predictionDir + 'pwms/'
	logodir = predictionDir + 'logos/'
	makeDir(predictionDir)
	makeDir(pwmdir)
	makeDir(logodir)
	
	# Get the directory containing the 3 position
	# experimental PWMs for this finger and set up
	# the output file for writing results
	expDir = '../data/revExp/F2_GAG/pwms3/'
	
	# Write header to results file
	fout = open(predictionDir + 'compare.txt', 'w')
	fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
	           %('num', 'targ','prot','canonprot','c1pcc','c2pcc',
	             'c3pcc', 'c1pcc.ic', 'c2pcc.ic', 'c3pcc.ic', 'p.c1cons',
	             'p.c2cons', 'p.c3cons', 'e.c1cons', 'e.c2cons', 'e.c3cons',
	             'pred.filt', 'finger', 'strin'))

	for fname in os.popen('ls ' + expDir):

		# Skip if this is the wrong stringency.
		fname = fname.strip()
		if strin == 'low' and re.match(r'(.)*_5mM.txt', fname) == None:
			continue
		
		# Get the info about this experiment
		sp_fname = fname.split('_')
		protNum = sp_fname[0]
		goal = sp_fname[1]
		prot = sp_fname[2].split('.')[0]
		canonProt = prot[0] + prot[2] + prot[3] + prot[6]
		label = '_'.join([str(protNum), goal, prot, strin])

		# Make the prediction for this protein
		nucMat = predict_matrix_bls([prot])

		# Make the nuc freq mat file and the logo
		makeNucMatFile(pwmdir, label, nucMat)
		logoIn = pwmdir + label + '.txt'
		logoOut = logodir + label + '.pdf'
		makeLogo(logoIn, logoOut, alpha = 'dna', 
		         colScheme = 'classic',
		         annot = "'5,M,3'",
		         xlab = '_'.join([goal,prot]))
		
		# Compare this pwm to the reverse experiment pwm
		expMat = pwmfile2matrix(expDir + fname)
		colPcc, colPcc_ic = comparePCC(nucMat, expMat)
		predCons = getConsensus(nucMat)
		#print predCons
		expCons = getConsensus(expMat)

		# Write the comparison results to file
		fout.write("%s\t%s\t%s\t%s\t" %(protNum, goal, prot, canonProt))
		fout.write("%.4f\t"*6 %(colPcc[0], colPcc[1], colPcc[2], \
		           				colPcc_ic[0], colPcc_ic[1], colPcc_ic[2]))
		fout.write("%s\t"*6 %(predCons[0], predCons[1], predCons[2],
		           			  expCons[0], expCons[1], expCons[2]))
		#print pred+'.'+filt
		fout.write("%s\t%s\t%s\n" %('bls02', 'F2', 'low'))

	fout.close()

def runMarcusDataAnalysis(style, decomp, weight_mat, order_mat, 
                          trainFing, trainStrin, useExact = True):
	# Runs the analysis of the for comparing the Marcus
	# experimental pwms using the given parameters

	# See how well BLS02 predicts the reverse experiments
	if style == 'bls02':
		runBLS02Analysis("../data/bls02RevExp/F2/low/","low")
		return

	# PErform the lookup or nn strategy on various datasets
	testFings = ['F2']
	testStrins = ['low']
	filts = ['cut10bc_0_5', 'cut3bc_0_5', 'cut10bc_0', 'cut3bc_025', 'cut10bc_025']
	filtsLabs = ['c10_0_5', 'c3_0_5', 'c10_0', 'c3_025', 'c10_025']
	for f in testFings:
		for s in testStrins:
			for i, filt in enumerate(filts):
				
				if re.match(r'top[0-9][0-9]', style) != None:
					topk = eval( style[(len(style) - 2):] )
				else:
					topk = None
				setWeightMatrices(weight_mat, order_mat)
				label, inDirPref = getLabels(style, decomp, weight_mat, useExact)
				decompDict = getDecompDict(style, decomp)

				inDir = '../data/b1hData/newDatabase/6varpos/' \
					+ trainFing + '/' + trainStrin + '/' + 'protein_seq_' + filt + '/'

				outDir = inDirPref + trainFing + '/' + trainStrin + \
					'/' + filt + '/'

				# Get the dictionary of binding frequencies
				canonical = True
				varpos = 6
				canInd = getPosIndex(varpos, canonical)
				freqDict = computeFreqDict(inDir, canInd)

				# Perform the actual prediction against the Marcus pwms
				if style == 'lookonly':
					lookupMarcusPWMs(inDir, outDir, freqDict, f, s, filtsLabs[i],
				    	             label, useNN = False,
				        	         skipExact = False, decompose = None)
				elif style == 'nnonly':
					lookupMarcusPWMs(inDir, outDir, freqDict, f, s, filtsLabs[i],
				                 	label, useNN = True, skipExact = True,
				                 	decompose = decompDict, topk = None)
				elif re.match(r'top[0-9][0-9]', style) != None and useExact:
					#print label
					lookupMarcusPWMs(inDir, outDir, freqDict, f, s, filtsLabs[i],
				                 	label, useNN = True, skipExact = False,
				                 	decompose = decompDict, topk = topk)	
				elif re.match(r'top[0-9][0-9]', style) != None:
					#print label
					lookupMarcusPWMs(inDir, outDir, freqDict, f, s, filtsLabs[i],
				                 	label, useNN = True, skipExact = True,
				                 	decompose = decompDict, topk = topk)	


def norm_matrix(matrix): #normalize each column to 1 -- Written by Anton Persikov
    columns = len(matrix) #get number of columns
    rows = len(matrix[0]) #get number of rows
    matrix_n = np.zeros((columns,rows), float)
    for i in range(columns):
        #read column
        summa = 0
        for j in range(rows):
            summa += matrix[i][j]
        #normalize column
        for j in range(rows):
            matrix_n[i][j] = float(matrix[i][j]) / summa
    return matrix_n

def predict_matrix_bls(protein): #predict matrix by BLS02 -- written by Anton Persikov
    # ref = Panayiotis V. Benos and Alan S. Lapedes and Gary D. Stormo, JMB, 2002, v.323(4): 701-727.
    # energy matrix was taken from table 2
    positions = ['03','04','02','01'] #positions in original order 
    ma_bls = {
        'A':[[ 1.29, 0.92, 1.72, 6.42], [ 0.25,-0.87, 1.32, 0.58], [ 0.54, 2.12, 2.27, -1.14], [-0.33, 1.38,-0.13, 0.72]],
        'C':[[ 5.21, 4.71, 6.29, 4.51], [ 3.86,-0.93, 6.92, 4.94], [ 4.16, 1.89, 6.27, -0.08], [ 4.68, 5.07, 5.11, 5.39]],
        'D':[[ 1.93,-2.35, 0.33,-0.68], [ 2.27, 0.42,-0.44,-0.29], [ 5.82,-1.32, 7.20,  6.32], [ 1.75, 0.64, 1.17,-0.21]],
        'E':[[ 0.34,-1.06,-0.05,-0.98], [ 6.34,-1.72, 1.88, 0.29], [ 1.30, 0.43, 1.98,  1.16], [ 0.09, 1.68, 1.26,-0.08]],
        'F':[[ 5.21, 4.71, 1.33, 4.51], [ 3.86, 5.55, 1.36, 4.94], [ 4.16, 6.92, 6.27,  5.13], [ 4.68, 0.42, 5.11, 5.39]],
        'G':[[ 2.45, 0.44, 1.07, 0.47], [ 0.23,-1.58, 0.39,-0.55], [ 0.27, 1.43, 2.27, -0.51], [ 2.62, 2.62, 1.61, 1.54]],
        'H':[[ 1.41,-0.10, 1.42,-0.85], [-0.38,-0.99, 0.20,-1.78], [ 0.28, 2.64,-1.31,  6.72], [ 1.73, 6.10,-0.08, 0.33]],
        'I':[[ 1.55, 1.58, 1.16, 1.77], [ 1.05, 5.55, 8.45, 5.90], [ 5.68, 2.58, 7.23,  1.23], [ 1.81, 1.77, 5.78, 6.46]],
        'K':[[ 0.33,-0.37,-1.04,-1.14], [ 6.59,-0.94, 2.70, 6.47], [ 5.54, 7.47, 0.77,  1.04], [ 0.87, 6.99,-1.01,-1.38]],
        'L':[[ 1.18, 0.19, 7.89,-0.46], [ 7.05, 6.63, 3.50, 1.14], [ 6.93, 1.54, 8.16,  0.89], [ 1.09, 2.44, 0.79, 7.94]],
        'M':[[ 6.52, 4.99, 6.55,-1.55], [ 6.33, 5.55, 2.41, 0.28], [ 0.26, 7.32, 2.06, -0.29], [ 6.91, 5.71, 5.48, 6.27]],
        'N':[[ 0.04, 0.73, 0.06,-1.33], [ 1.14,-1.72, 0.75, 1.21], [-3.73, 0.86,-0.005,-0.80], [-0.81, 0.48,-1.00, 1.75]],
        'P':[[ 7.33, 5.89, 7.36, 0.75], [ 7.01, 0.17, 2.27, 6.58], [ 6.20, 8.16, 2.85,  1.73], [ 7.70, 2.18, 0.21, 0.74]],
        'Q':[[-2.50, 6.32, 0.72,-0.52], [ 0.20,-1.32, 1.78, 0.30], [-0.26, 1.34, 0.77,  0.35], [-0.92, 6.10,-0.08, 0.33]],
        'R':[[ 1.01, 3.26,-1.36, 0.46], [ 1.66,-1.05, 2.11, 7.37], [ 6.59, 3.54, 2.56,  2.14], [ 0.38, 0.55,-3.33, 1.29]],
        'S':[[ 1.06, 0.73, 0.61,-0.57], [ 2.29,-1.45, 0.28,-0.91], [ 6.06, 1.15, 1.69, -1.65], [-0.21, 0.03,-0.40,-0.39]],
        'T':[[ 0.27, 1.72, 0.59,-1.33], [-0.18,-1.52, 1.29, 1.23], [ 7.50,-0.04, 3.26,  0.26], [-1.19, 0.20,-0.63,-0.46]],
        'V':[[ 7.59, 6.70, 1.57, 0.92], [ 7.01, 0.17, 1.76, 0.57], [ 8.01, 0.79, 8.52,  0.80], [ 1.13, 1.98, 0.46, 0.78]],
        'W':[[ 2.21, 4.71, 6.29, 4.51], [ 3.86,-0.52, 0.67,-0.71], [ 4.16, 6.92, 6.27, -0.08], [ 4.68, 5.07, 5.11,-0.27]],
        'Y':[[ 5.21, 4.71,-0.05, 4.51], [ 4.95,-1.25, 1.43, 5.31], [ 4.16, 6.92, 1.27,  5.13], [-1.37, 5.85, 5.23, 5.61]]}
    aa_pos = {'01':6, '02':3, '03':0, '04':2} #canonical amino acid positions
    Nzf = len(protein)
    #pwm = np.zeros((Nzf*3+1,4), float) #PWM for the DNA length Nzf*3+1
    pwm = np.zeros((Nzf*3,4), float) #PWM for the DNA length Nzf*3
    for z in range(Nzf):
        nuc1 =((Nzf-1)-z)*3 #first base in corresponding DNA sub-sequence
        for i in range(4): #realtive nucleotide position
            pos = '0'+repr(i+1) #contacting position
            a = protein[z][aa_pos[pos]] #contacting amino acid
            for j in range(4):
                prob = np.exp(-ma_bls[a][positions.index(pos)][j]) #probability for current position
                if pos == '04': #junction position
                	# Skip the 4th base position since we only have 3 positions
                	continue
                	"""
                    if z != 0: #there are probabilities from previous ZF
                        pwm[nuc1+i][j] = pwm[nuc1+i][j] * prob #use combined probability
                    else:
                        pwm[nuc1+i][j] = prob #no complementarity for BLS!
                    """
                else:
                    pwm[nuc1+i][j] = prob
    return norm_matrix(pwm)

def main():

	#styles = ['bls02']
	#styles = ['nnonly']
	#styles = ['top20', 'top25', 'top30', 'top35', 'top40']
	styles = ['top10', 'top15']
	#styles = ['top30']
	weight_mats = [None, 'PAM30']
	useExact = True
	#weight_mats = [None]#['PAM30']
	#decomps = ['singles', 'doubles', 'triples']
	decomps = ['single']
	order_mat = 'PAM30'
	trainFing = "F2"
	trainStrin = "low"

	# Run the "topk" neighbors analysis
	for style in styles:
		for weight_mat in weight_mats:
			for decomp in decomps:
				print "Running:\t%s\t%s\t%s\t%s" %(style, decomp, weight_mat, 
				                                   order_mat)
				runMarcusDataAnalysis(style, decomp, weight_mat, order_mat, 
                          		  	  trainFing, trainStrin, useExact)

if __name__ == '__main__':
	main()