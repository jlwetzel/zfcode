# Script for manipulating parsing through the B1H database
# starting at the level of variable protein regions in the 
# form of codons.

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

def parseMasterFile(path):
	# Returns a dictionary:
	# (seqrun, barcode) -> (target, finger, stringency)
	mfile = open(path, 'r')
	bcodes = {}
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
					bcodes[(seqrun, barcode)] = (target, finger, 
					                             stringency)
	mfile.close()
	return bcodes

def main():
	mFilePath = '../data/b1hdata/database/all_selection_data_final.txt'
	bcodes = parseMasterFile(mFilePath)
	for key in sorted(bcodes.keys()):
		print str(key) + ":" + str(bcodes[key])

if __name__ == '__main__':
	main()