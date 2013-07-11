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

import os
import math
import numpy as np

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

def outputSeqStats(path, lineProc, nonNNS, totalCount):
    # Output statsistics about number of sequences processed,
    # number of sequencees thrown away, and total counts.

    fout = open(path, 'w')
    fout.write('Unique sequences containing only NNS codons: %d\n'
                %(lineProc - nonNNS))
    fout.write('# seqs that contained nonNNS codons: %d\n'
                %(nonNNS))
    fout.write('Total count accross all kept (NNS) seqs: %d\n'
                %totalCount)
    fout.close()


def getStatsAndFilterNNS(path):
    # Processes each all_nuc_seq.txt file in the directory 
    # given by path.  A new directory is created which 
    # which contains one file for each file in path, 
    # with sequences contianing nonNNS codons removed
    # and counts taken log2(count + 2) then normalized per file.

    nucs = ['A', 'C', 'T', 'G']

    # Make a directory for nucleotide and codon stats
    newDir = '/'.join(path.split('/')[:-2]) + '/all_nuc_seq_NNSnorm/'
    try:
        os.mkdir(newDir)
    except OSError:
        pass
    try:
        os.mkdir(newDir + "statistics/")
    except OSError:
        pass

    # For each file beginning with an capital nucleotide letter
    handle = os.popen("ls " + path, 'r')
    k = 0
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
        totalCount = 0
        
        for line in fin:
            sp_line = line.strip().split()
            totalCount += int(sp_line[1])
            seq, count = sp_line[0], \
                math.log(int(sp_line[1]) + 2,2)

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
            updateSeqCounts(seq, count, fileSeqCounts)
            updateCodonCounts(seq, count, fileCodonCounts)
            updateNucCounts(seq, count, fileNucCounts)

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
        fin.close()
        

def main():
    fings = ['F1', 'F2', 'F3']
    strins = ['low', 'high']
    for fing in fings:
        for strin in strins:
            path = '../data/b1hData/newDatabase/6varpos/' + \
                fing + '/' + strin + '/all_nuc_seq/'
            getStatsAndFilterNNS(path)

if __name__ == '__main__':
    main()



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
#     compute the JS divergece these codon frequencies from 
#     the observed codon usage in E. coli given in the paper
#     by Malloy, Stewart, and Taylor, 1996.  (We first normalize
#     these observed frequencies to get rid of non-NNS codons.)