# A module for creating pwm files.

def initPosCounts(npos, type):
	# Initialize the dictionay to all 0 counts
	nucs = ['A', 'C', 'G', 'T']           
	aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

	for p in range(npos):
		if type == 'protein':
			for a in aminos:
				posCounts[(p,a)] = 0
		elif type == 'dna':
			for n in nucs:
				posCounts[(p,n)] = 0

def makePwm(src, dst, npos, type):
	# src is a tab delimited txt file where the first column 
	# is a list of proteins or dnas (each of same length)
	# and the second is a list of counts or frequencies.
	# Outputs a pwm file in the transfac format to dst.
	#
	# npos is the number of positions per sequence
	# type is either 'protein' or 'dna'

	posCounts = {}  # (pos, letter) -> count
	initPosCounts(npos, type)
	letters = sorted([i[1] for i in posCounts.keys()])
	print letters

	# Fill the dictionary
	fin = open(src, 'r')
	for seq in fin:
		pass

def main():
	src = '../data/b1hData/newDatabase/6varpos/' + \
		'F2/high/protein_cut10_entr04/AAA.txt'
	makePwm('tmp.txt', '', 6, "protein")

if __name__ == '__main__':
	main()