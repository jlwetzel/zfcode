# A script for converting the *all_nuc_seq.txt files from 
# codons of variable protein positions to their to new files
# that contain the proteins encoded by those codons.
# Additionally, the script outputs a variety of statistics 
# about each *all_nuc_seq.txt file (e.g. codon usage).
#
# Outputs new files, with only one file per 3-mer target,
# stringency, and finger combination with information 
# about the number of ways each protein was coded and the 
# JS divergence for the distrubution of codons from the 
# expected distribution of codons for E. coli (if using 
# only NNS codons).  
#
#
# Processing conventions:
#
# 1.  Any sequence containing any non-NNS codons will be 
#     assumed a sequencing error and thrown out.
# 2.  We take the log2 of the counts (+2) for each codon combination
#     assuming count to be proportional to bacterial growth rate.
# 3.  The log-counts for each codon combination in a given 
#     *all_nuc_seq.txt file are converted to frequencies, as we 
#     are really only concerned with relative affinity of different 
#     proteins for a given 3-mer target.
# 4.  When combining multiple files from different sequencing
#     runs for the same 3-mer target, we take an arithmetic 
#     mean for the count-derived freqencies of each codon 
#     combination.
# 5.  Based on the codon frequencies derived in step 4, we 
#     compute the JS divergece of these codon frequencies from 
#     the observed codon usage in E. coli given in the paper
#     by Malloy, Stewart, and Taylor, 1996.  (We first normalize
#     these observed frequencies to get rid of non-NNS codons.)

import os
import math
import numpy as np
import itertools
import re

###################################################
###################################################
#    Global Constants
#
###################################################
###################################################

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


###################################################
###################################################
#    Helper Functions
#
###################################################
###################################################

def updateCodonCounts(seq, count, codonCounts):
    # Updates a dictionary of form: (codonPos/codon) -> count
    # based on the sequence and count
    
    i = 0
    codonNum = 0
    while i < len(seq):
        if i%3 == 2:
            codon = seq[i-2:i+1]
            if codonCounts.has_key((codonNum,codon)):
                codonCounts[codonNum,codon] += count
            else:
                codonCounts[codonNum,codon] = count
            codonNum += 1
        i += 1
    
def updateNucCounts(seq, count, nucCounts):
    # Updates a dictionary of form: (nucPos/codon)-> count
    # based on the sequence and count
    
    i = 0
    while i < len(seq):
        if nucCounts.has_key((i, seq[i])):
            nucCounts[i, seq[i]] += count
        else:
            nucCounts[i, seq[i]] = count
        i += 1

def updateSeqCounts(seq, count, seqCounts):
    # Updates a dictionary of form: sequence -> count
    # based on the sequence and count
    
    if seqCounts.has_key(seq):
        seqCounts[seq] += count
    else:
        seqCounts[seq] = count

def outputCodonStats(path, codonCounts, numNucs):
    # Outputs a table of form (codonPos/codon) -> count
    # where 1 row for each codon and one column for each 
    # codonPos.  Table is output to file 'path'.
    # numNucs is the number of nucleotides per sequence 
    # for this particular dataset.
    
    fout = open(path, 'w')
    numCodons = numNucs/3

    # Normalize each column of the data
    tempList = []
    for i in range(numCodons):
        keys = [(pos,codon) for (pos, codon) in codonCounts.keys() 
                if int(pos) == i]
        vals = np.array([codonCounts[k] for k in keys])
        vals = vals/vals.sum()
        for i in range(len(keys)):
            codonCounts[keys[i]] = vals[i]

    # print the table header
    header = "Cod\t"
    for i in range(numCodons - 1):
        header += 'c' + str(i) + '\t'
    header += 'c' + str(i + 1) + '\n'
    fout.write(header)

    # one line for each NNS codon
    for n1 in ['A', 'C', 'G', 'T']:
        for n2 in ['A', 'C', 'G', 'T']:
            for n3 in ['C', 'G']:
                codon = n1 + n2 + n3
                line = codon + '\t'
                for i in range(numCodons - 1):
                    try:
                        line += str(codonCounts[i, codon]) + '\t'
                    except KeyError:
                        line += '0\t'
                try:
                    line += str(codonCounts[i + 1, codon]) + '\n'
                except KeyError:
                    line += '0\n'
                fout.write(line)
    fout.close()

def outputNucStats(path, nucCounts, numNucs):
    # Outputs a table of the form (nucPos/nuc) -> count
    # where 1 row for each nuc and and 1 column for
    # each nuc position. 
    # numNucs is the number of nucleotides per sequence 
    # for this particular dataset. 
    
    fout = open(path, 'w')

    # Normalize each column of the data
    tempList = []
    for i in range(numNucs):
        keys = [(pos,nuc) for (pos, nuc) in nucCounts.keys() \
                if int(pos) == i]
        vals = np.array([nucCounts[k] for k in keys])
        vals = vals/vals.sum()
        for i in range(len(keys)):
            nucCounts[keys[i]] = vals[i]

    # print the table header
    header = "Nuc\t"
    for i in range(numNucs - 1):
        header += 'n' + str(i) + '\t'
    header += 'n' + str(i + 1) + '\n'
    fout.write(header)

    # One line for each nuc
    for n in ['A', 'C', 'G', 'T']:
        line = n + '\t'
        for i in range(numNucs - 1):
            try:
                line += str(nucCounts[i, n]) + '\t'
            except KeyError:
                line += '0\t'
        try:
            line += str(nucCounts[i + 1, n]) + '\n'
        except KeyError:
                line += '0\n'
        fout.write(line)
    fout.close()

def outputSeqStats(path, lineProc, nonNNS):
    # Output statsistics about number of sequences processed,
    # number of sequencees thrown away, and total counts.

    fout = open(path, 'w')
    fout.write('Unique sequences containing only NNS codons: %d\n'
                %(lineProc - nonNNS))
    fout.write('# seqs that contained nonNNS codons: %d\n'
                %(nonNNS))
    fout.close()

###################################################
###################################################
#    Step 1 Functions --
#    Filtering and stats collection.
###################################################
###################################################

def outputNNSNormFile(path, seqCounts):
    # Normalizes counts for all nuc seqs in the file to 
    # 1 and outputs them to the file named 'path'
    
    fout = open(path, 'w')
    sortedSeqs = sorted(seqCounts.items(), 
                        key=lambda x: (x[1],x[0]), reverse=True)
    seqs = [i[0] for i in sortedSeqs]
    counts = np.array([i[1] for i in sortedSeqs])
    counts = counts/counts.sum()
    for i in range(len(seqs)):
        fout.write(seqs[i] + '\t' + str(counts[i]) + '\n')
    fout.close()



def getStatsAndFilterNNS(path):
    # The first level of processing.
    #
    # Processes each all_nuc_seq.txt file in the directory 
    # given by path.  A new directory is created which 
    # which contains one file for each file in path, 
    # with sequences containing nonNNS codons removed
    # and counts taken log2(count + 2) then normalized per file.

    nucs = ['A', 'C', 'T', 'G']

    # Make directories for nucleotides sequences with 
    # normalized probs
    newDir = '/'.join(path.split('/')[:-2]) + '/all_nuc_seq_NNSnorm/'
    try:
        os.mkdir(newDir)
    except OSError:
        pass
    # Make a directory for nucleotide and codons stats. 
    try:
        os.mkdir(newDir + "statistics/")
    except OSError:
        pass

    
    handle = os.popen("ls " + path, 'r')
    for line in handle:

        # This is a subdirectory
        if line[0] not in nucs:
            continue

        # Gather info from the filename
        fname = line.strip()
        sp_fname = fname.split('_')
        targ, seqrun, bcode = \
            sp_fname[0], sp_fname[1], sp_fname[2]

        fileCodonCounts, fileSeqCounts, fileNucCounts = \
            {}, {}, {}

        print "Processing: %s" %(path + fname)
        # open the file and begin parse
        fin = open(path + fname, 'r')
        lineProc = 0
        nonNNS = 0
        
        for line in fin:
            sp_line = line.strip().split()
            count = int(sp_line[1])
            seq, logcount = sp_line[0], \
                math.log(count + 2,2)

            # Does seq contain only NNS codons?
            lineProc += 1
            keep = True
            for i in range(len(seq)):
                if i%3 == 2 and seq[i] in ['A', 'T']:
                    keep = False
                    break
            if not keep:
                nonNNS += 1
                continue

            # Update statistics for this file
            updateSeqCounts(seq, logcount, fileSeqCounts)
            updateCodonCounts(seq, logcount, fileCodonCounts)
            updateNucCounts(seq, logcount, fileNucCounts)

        # Output the new files
        newfname = fname.split('.')[0] + '_NNS_norm.txt'
        outputNNSNormFile(newDir + newfname, fileSeqCounts)
        newfname = fname.split('.')[0] + '_codonStats.txt'
        outputCodonStats(newDir + 'statistics/' + newfname,
                         fileCodonCounts, len(seq))
        newfname = fname.split('.')[0] + '_nucStats.txt'
        outputNucStats(newDir + 'statistics/' + newfname,
                       fileNucCounts, len(seq))
        newfname = fname.split('.')[0] + '_seqStats.txt'
        outputSeqStats(newDir + 'statistics/' + newfname,
                       lineProc, nonNNS)
        fin.close()

###################################################
###################################################
#    Step 2 Functions --
#    Combining frequencies of multiple experiments
###################################################
###################################################

def outputCombinedNucFile(path, targ, seqDict, numTargFiles):
    # Output a single file with the average frequency
    # for each unique nuc sequence of seqDict.

    # Get the mean for the list of frequencies for 
    # each sequence
    for seq in seqDict.keys():
        seqDict[seq] = \
            np.array(seqDict[seq]).sum()/float(numTargFiles)

    # Sort the sequences by decreasing value
    sortedSeqs = sorted(seqDict.items(), 
                        key=lambda x: (x[1],x[0]), reverse=True)
    seqs = [i[0] for i in sortedSeqs]
    freqs = [i[1] for i in sortedSeqs]

    # Output combined seqs and freqs to the new file
    fout = open(path + targ + '_combined_nuc_seq.txt', 'w')
    for i in range(len(seqs)):
        fout.write(seqs[i] + '\t' + str(freqs[i]) + '\n')
    fout.close()
    print 'Output results to %s' \
        %(path + targ + '_combined_nuc_seq.txt')

def combineExperiments(path):
    # Combines the multiple experiments for into one file
    # per 3-mer by taking the arithmetic mean of frequencies
    # for each sequence of nucleotides across the experiments.

    nucs = ['A', 'C', 'T', 'G']

    # Make a directory for nucleotide and codon stats
    newDir = '/'.join(path.split('/')[:-2]) + '/combined_nuc_seq/'
    try:
        os.mkdir(newDir)
    except OSError:
        pass
    
    targ = '' # current 3-mer target
    numTargFiles = 1 # number of files for the current 3mer targ

    handle = os.popen("ls " + path, 'r')
    k = 0
    for line in handle:

        print line
        # This is a subdirectory
        if line[0] not in nucs:
            continue

        # Gather info from the filename
        prevTarg = targ       # 3-mer target for the last file
        fname = line.strip()
        sp_fname = fname.split('_')
        targ, seqrun, bcode = \
            sp_fname[0], sp_fname[1], sp_fname[2]

        # Create a dict for mapping sequences to list of 
        # frquencies for that sequence across the 
        # experiments.
        
        if targ != prevTarg:
            
            if k > 0:
                outputCombinedNucFile(newDir, prevTarg, 
                                      seqDict, numTargFiles)
            seqDict = {}
            numTargFiles = 1
        else:
            numTargFiles += 1

        # Open the current file and process
        fin = open(path + fname, 'r')

        print "Processing: %s" %(path + fname)
        for line in fin:
            sp_line = line.strip().split()
            seq = sp_line[0]
            freq = eval(sp_line[1])

            if seqDict.has_key(seq):
                seqDict[seq].append(freq)
            else:
                seqDict[seq] = [freq]

        fin.close()
        k += 1
    
    # For the last target encountered.
    outputCombinedNucFile(newDir, targ, seqDict, numTargFiles)

###################################################
###################################################
#    Step 3 Functions --
#    Convert to proteins and get JS divergence
#    from background E coli codon usage
###################################################
###################################################

def getCodonBiasNNS():
    # Returns a dictionary of dictionaries for codon bias
    # in e. coli arranged by the amino acid that the codon
    # codes for (see below).  Non-NNS codons are removed, then
    # then distributions are renormalized.

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

    nnsBias = {}
    for amino in codonBias:
        nnsBias[amino] = {}
        codonDict = codonBias[amino]
        items = [i for i in codonDict.items() \
                if i[0][2] in ['C', 'G']]
        codons = [i[0] for i in items]
        freqs = np.array([i[1] for i in items])
        freqs = freqs/freqs.sum()
        for i in range(len(codons)):
            nnsBias[amino][codons[i]] = freqs[i]

    return nnsBias

def klDiv(f, g):
    # Computes the KL divergence for discrete 
    # distributions f and g, which are lists of 
    # frequencies for f and g (in the same order)
    # g must be non-zero at all points of the PMF
    # for KL to be defined.
    #
    # f[i]lg(f[i]/g[i]) is taken to be 0 if f[i] = 0.

    kl = 0.0
    for i in range(len(f)):
        if f[i] == 0:
            continue
        else:
            kl += f[i]*math.log(f[i]/g[i], 2)

    return kl

def jsDiv(f, g):
    # Computes the Jensen Shannon divergence of 
    # f and g, where f and g are lists for discrete 
    # PMFs (listed in the same order)

    m = 0.5*(f + g)
    jsd = 0.5*(klDiv(f, m) + klDiv(g, m))
    return jsd

def getProtJSDiv(prot, codonDict, nnsBiasEcoli):
    # Returns the JS diveregnce for the frequency of
    # different codon combinations for a protein 
    # from the seltection vs the expected background
    # freqeuncies for those coding combinations 
    # based on E. coli codon usage
    #
    # prot is the protein
    # codon dict is a dictionary mapping codon 
    #  combinations to frequencies in the selection
    # see getCodonBiasNNS for def of nnsBiasEcoli

    # Make a list of lists of codons for each position
    # of the protein.
    
    protCodons = []
    for i, amino in enumerate(prot):
        protCodons.append([])
        for codon in nnsBiasEcoli[amino].keys():
            protCodons[i].append(codon)

    # Find expected background frquency for each codon
    # combination assuming independence of codons 
    comboList = list(itertools.product(*protCodons))
    comboBGProbs = []
    for combo in comboList:
        prob = 1
        for i, codon in enumerate(combo):
            prob *= nnsBiasEcoli[prot[i]][codon]
        comboBGProbs.append(prob)
    comboBGProbs = np.array(comboBGProbs)

    # Get the list of frequencies for codon combinations
    # found in the selection, order the same is the
    # combination bg probs
    comboList = [''.join(i) for i in comboList]
    obsComboProbs = []
    
    for combo in comboList:
        if codonDict.has_key(combo):
            obsComboProbs.append(codonDict[combo])
        else:
            obsComboProbs.append(0.0)

    f = np.array(obsComboProbs)
    g = np.array(comboBGProbs)
    f = f/f.sum()
    g = g/g.sum()

    # Return the (#possible combos, JSD) of the distributions
    return len(comboList), jsDiv(f, g)

def convertToProteins(path):
    # Converts a directory with files of nucleotide
    # sequences and frequencies to the corresponding
    # combined set of proteins and frequencies.  Also JS
    # divergence is calculated for the ways the protein
    # is coded in the selections compared to E coli 
    # bg codon usage.  This info is output into a 
    # new directory with corresponding files.

    nucs = ['A', 'C', 'T', 'G']
    nnsBiasEcoli = getCodonBiasNNS()

    # Make directories for nucleotides sequences with 
    # normalized probs
    newDir = '/'.join(path.split('/')[:-2]) + '/protein_seqs_JSD/'
    try:
        os.mkdir(newDir)
    except OSError:
        pass
    # Make a directory for nucleotide and codons stats. 
    try:
        os.mkdir(path + "statistics/")
    except OSError:
        pass

    
    handle = os.popen("ls " + path, 'r')
    for line in handle:

        # This is a subdirectory
        if line[0] not in nucs:
            continue

        # Gather info from the filename
        fname = line.strip()
        sp_fname = fname.split('_')
        targ, seqrun, bcode = \
            sp_fname[0], sp_fname[1], sp_fname[2]

        # Create dictionary of dictionaries mapping prot seq
        # to dict of observed codon sequences which are 
        # in turn mapped to the observed frequencies of 
        # these codon sequences.  Also make a dictionary
        # mapping (position/codon) -> codon frequency
        protDict = {} 
        codonPosFreq = {}
        nucPosFreq = {}

        # Open the current file and process
        fin = open(path + fname, 'r')
        print "Processing: %s" %(path + fname)
        for line in fin:
            sp_line = line.strip().split()
            nucseq, freq = sp_line[0], eval(sp_line[1])

            #print nucseq, freq
            # Build the protein sequence
            prot = ''
            codonpos = 0
            for i in range(len(nucseq)):
                if i%3 == 2:
                    # Grab the codon and update codonPosFreq
                    codon = nucseq[i-2:i+1]
                    prot += codon2amino[codon]
                    if codonPosFreq.has_key((codonpos, codon)):
                        codonPosFreq[codonpos,codon] += freq
                    else:
                        codonPosFreq[codonpos,codon] = freq
                if nucPosFreq.has_key((i, nucseq[i])):
                    nucPosFreq[i, nucseq[i]] += freq
                else:
                    nucPosFreq[i, nucseq[i]] = freq
                    codonpos += 1

            # Update the protein/nucseq -> freq dictionary
            if protDict.has_key(prot):
                protDict[prot][nucseq] = freq
            else:
                protDict[prot] = {}
                protDict[prot][nucseq] = freq
            
        fin.close()

        print "Getting JS divergence for proteins."
        # Create list of (freq, prot, #observed combos,
        # possible combos, jsd) tuples
        tList = []
        #print len(protDict.keys())
        for prot in protDict.keys():
            #print len(protDict[prot])
            numPos, jsd = getProtJSDiv(prot, protDict[prot],
                                       nnsBiasEcoli)
            totProtFreq = 0.0
            for s in protDict[prot]:
                totProtFreq += protDict[prot][s]
            #print totProtFreq
            tList.append((totProtFreq, prot, 
                          len(protDict[prot]), 
                          numPos, jsd))
        
        #print protDict['WRSWLA']  
        print "Sorting by frequency."  

        # Sort the tuple list by frequency and output it
        tList = sorted(tList, reverse = True)

        print "Outputting to files."
        fout = open(newDir + targ + '_protein_seqs_JSD.txt', 'w')
        for x in tList:
            fout.write(x[1] + '\t' + str(x[0]) + '\t' \
                       + str(x[2]) + '\t' + str(x[3]) + '\t' \
                       + str(x[4]) + '\n')

        # Output the codon usage statistics
        newfname = fname.split('.')[0] + '_codonStats.txt'
        outputCodonStats(path + './statistics/' + newfname,
                         codonPosFreq, len(nucseq))
        # Output the nucleotide usage statistics
        newfname = fname.split('.')[0] + '_nucStats.txt'
        outputNucStats(path + './statistics/' + newfname,
                         nucPosFreq, len(nucseq))


###################################################
###################################################
#    MAIN
#
###################################################
###################################################

def main():

    # For 6 variable positions
    fings = ['F3', 'F2', 'F1']
    strins = ['low', 'high']
    for fing in fings:
        for strin in strins:
            # Step 1 -- Filter and gather stats.
            path = '../data/b1hData/newDatabase/6varpos/' + \
                fing + '/' + strin + '/all_nuc_seq/'
            getStatsAndFilterNNS(path)

            # Step 2 -- Combine multiple experiment frequencies
            path = '../data/b1hData/newDatabase/6varpos/' + \
                fing + '/' + strin + '/all_nuc_seq_NNSnorm/'
            combineExperiments(path)

            # Step 3 -- Convert to proteins and compute
                        JS divergence for each protein
            path = '../data/b1hData/newDatabase/6varpos/' + \
                    fing + '/' + strin + '/combined_nuc_seq/'
            convertToProteins(path)
            
    # For 5 variable positions
    fings = ['F2', 'F3']
    strins = ['low', 'high']
    for fing in fings:
        for strin in strins:
            
            # Step 1 -- Filter and gather stats.
            path = '../data/b1hData/newDatabase/5varpos/' + \
                fing + '/' + strin + '/all_nuc_seq/'
            getStatsAndFilterNNS(path)

            # Step 2 -- Combine multiple experiment frequencies
            path = '../data/b1hData/newDatabase/5varpos/' + \
                fing + '/' + strin + '/all_nuc_seq_NNSnorm/'
            combineExperiments(path)

            # Step 3 -- Convert to proteins and compute
                        JS divergence for each protein
            path = '../data/b1hData/newDatabase/5varpos/' + \
                    fing + '/' + strin + '/combined_nuc_seq/'
            convertToProteins(path)

if __name__ == '__main__':
    main()
