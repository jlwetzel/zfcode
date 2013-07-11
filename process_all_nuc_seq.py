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
# Processing conventions:
#
# 1.  Any sequence containing any non-NNS codons will be 
#     assumed a sequencing or contamination error and thrown out.
# 2.  We take the log2 of the counts for each codon combination
#     assuming count to be proportional to bacterial growth rate.
# 3.  If the total log2 count added across all codon combinations 
#     forcodon a given protein is less than 3 we throw away all 
#     those codon combinations, under the assumption that 
#     non-specific binding has occurred.
# 4.  After steps 2 and 3 are completed, the log-counts for 
#     each codon combination in a given *all_nuc_seq.txt file
#     are converted to frequencies, as we are really only 
#     concerned with relative affinity of different proteins 
#     for the corresponding 3-mer targets.
# 5.  When combining multiple files from different sequencing
#     runs for the same 3-mer target, we take an arithmetic 
#     mean for the count-derived freqencies of each codon 
#     combination.
# 6.  Based on the codon frequencies derived in step 6, we 
#     compute the JS divergece these codon frequencies from 
#     the observed codon usage in E. coli given in the paper
#     by Malloy, Stewart, and Taylor, 1996.  (We first normalize
#     these observed frequencies to get rid of non-NNS codons.)

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

