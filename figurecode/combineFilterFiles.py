fing = "F2"
strin = "low"

filtDirPref = '../../data/b1hData/newDatabase/6varpos/' + \
	fing + '/' + strin + '/'
filters = ['cut10bc_0_5', 'cut10bc_025', 'cut3bc_0_5', 'cut3bc_025']
filtDirs = [(filtDirPref + '_'.join(["protein", "seq", x])) for x in filters]
filtFiles = ['/'.join([x, 'all.csv']) for x in filtDirs]

fout = open('../../figures/aminoFreqs/F2_low_combinedFilters.txt', 'w')
for i in range(len(filtFiles)):
	fin = open(filtFiles[i], 'r')
	if i != 0:
		fin.readline()  # strip the header
	for j, l in enumerate(fin):
		l = l.strip().split('\t')
		if i == 0 and j == 0:
			l = ','.join(l + ["filter"]) + '\n'
		else:
			l = ','.join(l + [filters[i]]) + '\n'
		fout.write(l)
	fin.close()
fout.close()