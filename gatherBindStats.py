# Basic analysis of the diverity of the binding proteins.

import os
import sys
import math
import re
from itertools import product

nucs = ['A', 'C', 'G', 'T']
aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def getUniqueShilpaZFs(path, contacts):
	# Returns a list of ZFs for the organism of the specified
	# file (when using only given list of contact positions)
	# Assumes are in Shilpa's format and have domains listed 
	# from -1,1,2,3,4,5,6,7 positions.

	inFile = open(path, 'r')
	lines = [l.strip().split() for l in inFile if l[0] != '#']
	regions = [l[11] for l in lines if l[12] == "NOGAP"]
	prots = set()
	for r in regions:
		prot = ""
		for i in contacts:
			prot += r[i]
		prots.add(prot)

	return prots

def getProtDict(path, contacts):
	# Returns a dictionary mapping unique proteins 
	# (when using only given list of contact positions)
	# to a list of tuples (3-mer target bound, count).
	# path should be to a file the 'all.txt' format
	# and contatcs positions in file are -1,1,2,3,5,6
		
	inFile = open(path, 'r')
	inFile.readline()  # Skip the header line
	protDict = {}

	for line in inFile:
		sp_line = line.strip().split('\t')
		freq = eval(sp_line[2])
		targ = sp_line[0]
		prot = ""
		for c in contacts:
			prot += sp_line[1][c]
            
		if protDict.has_key(prot):
			if protDict[prot].has_key(targ):
				protDict[prot][targ] += freq
			else:
				protDict[prot][targ] = freq
		else:
			protDict[prot] = {}
			protDict[prot][targ] = freq

	inFile.close()
	return protDict

def getTargDict(path, contacts, NNpos):
    # Returns a dictionary mapping each 3-mer target 
    # to a set of proteins that are bound by that target
    # Assumes path is a directory name for a collection
    # of protein files with names like AAA.txt, AAC.txt
    # where the first column is a list of proteins bound.
    #
    # Contatcs is the list of indices to be considered 
    # when forming the set of proteins.
    # (i.e. if contatcs = [0,2,3,5] then any two proteins
    #  that have the same amino acids at these positions
    #  are considered to be identical)
    # 
    # NNpos is set to False by default, but if a list of 
    # indices is passed to NNpos, then the protein lists
    # with be augmented with nearest neighbors (i.e. proteins
    # that are off by only one amino acid with respect to 
    # the indices passed to NNpos).  These indicies should
    # be given with respect to the protein that will 
    # that would be left after converting according to 
    # the indices given by the contacts variable.

    targDict = {}
    for fname in os.popen('ls ' + path):
        fname = fname.strip()
        
        # Skip files that aren't labeled in the expected way
        if re.match(r'[ACGT]{3}.txt', fname) == None:
            continue

        # Set up the dictionary for this 3mer target
        targ = fname.split('.')[0]
        targDict[targ] = set()
        fin = open(path + fname, 'r')

        # Fill in the dicitonare for this 3-mer target
        for line in fin:
            sp_line = line.strip().split()
            origProt = sp_line[0]
            newProt = ''
            for i in contacts:
                newProt += origProt[i]
            targDict[targ].add(newProt)

            # Allow nearest neighbors if requested
            if NNpos == None:
                continue
            for i in NNpos:
                for a in aminos:
                    if newProt[i] == a:
                        continue
                    neighbor = newProt[:i] + aminos.index(a) + \
                        newProt[i+1:]
                    targDict[targ].add(neighbor)
        fin.close()

    return targDict


def getProtSet(path, contacts):
	# Returns the set of unique proteins for the 
	# contact positions given by contacts.
	# Path should be to a file of the 'all.txt' format.
	# Only includes those proteins which have a JSD 
	# of at most jsdCut from the expected E. coli 
	# codon distribution for that protein.
	
	inFile = open(path, 'r')
	inFile.readline()  # Skip the header line
	protSet = set()

	for line in inFile:
		sp_line = line.strip().split('\t')
		prot = ""
		for c in contacts:
			prot += sp_line[1][c]
		protSet.add(prot)

	inFile.close()
	return protSet

def printBindingSetStats(fings, strins, bindset, maxSize, allProts = None):
    
    # Get the sizes of the individual sets.
    numpos = math.log(maxSize, 20)
    print "Max Size for %d positions: %d" %(numpos, maxSize)
    print "#"*64
    for f in fings:
        for s in strins:
            setSize = len(bindset[f,s])
            print "Size of %s %s: %d" %(f, s, setSize) 
        
    # Check if high stringency is a subest of low stringency.
    print
    print "Is high a subset of low? (|High - Low|/|High|)"
    print "#"*64
    for f in fings:
        highMinLow = bindset[f,'high'] - bindset[f,'low']
        print "%s: %.3f" %(f, len(highMinLow)/ \
                           float(len(bindset[f,'high'])))
    
    # Get sizes for unions/intersections of different combinations
    # of fingers.
    print 
    print "Combining Fingers"
    print "#"*64
    for i in range(len(fings)):
        for j in range(i + 1, len(fings)):
            for s in strins:
            	inter = bindset[fings[i], s] & bindset[fings[j], s]
            	union = bindset[fings[i], s] | bindset[fings[j], s]
            	print "%s %s Intersection (%s): %d" \
            		%(fings[i], fings[j], s, len(inter))
            	print "%s %s Union (%s)       : %d" \
            		%(fings[i], fings[j], s, len(union))
            	print "%s %s Jaccard (%s)     : %.3f" \
            		%(fings[i], fings[j], s, len(inter)/float(len(union)))
        	print
    for s in strins:
    	inter = bindset[fings[0], s]
    	union = set()
    	for f in fings:
    		union = union | bindset[f,s]
    		inter = inter & bindset[f,s]
    	print "All Intersection (%s): %d" %(s, len(inter))
        print "All Union (%s): %d" %(s, len(union))

    # Get sizes for combining high and low stringencies.
    print
    print "Combining Stringencies"
    print "#"*64
    bothStrinsInter = {}
    bothStrinsUnion = {}
    for i in range(len(fings)):
        inter = bindset[fings[i], 'low'] & bindset[fings[i], 'high']
        bothStrinsInter[fings[i]] = inter
        union = bindset[fings[i], 'low'] | bindset[fings[i], 'high']
        bothStrinsUnion[fings[i]] = union
        print "All of %s Intersection: %d" \
                    %(fings[i], len(inter))
        print "All of %s Union       : %d" \
                    %(fings[i], len(union))
    
    # Unioning the intersections for both stringencies
    print 
    print "Unioning the intersections"
    unionInter = {}
    for i in range(len(fings)):
        for j in range(i+1, len(fings)):
            union = bothStrinsInter[fings[i]] | bothStrinsInter[fings[j]]
            unionInter[fings[i], fings[j]] = union
            print "All of %s and %s: %d" \
                %(fings[i], fings[j], len(union))
    if len(fings) == 3:
        union = bothStrinsInter['F1'] | bothStrinsInter['F2'] | bothStrinsInter['F3']
        unionInter['F1', 'F2', 'F3'] = union
        print "All of F1 F2 F3 : %d" %len(union)

    # Intersecting the union for both stringencies
    print 
    print "Intersecting the unions"
    interUnion = {}
    for i in range(len(fings)):
        for j in range(i+1, len(fings)):
            inter = bothStrinsUnion[fings[i]] & bothStrinsUnion[fings[j]]
            interUnion[fings[i], fings[j]] = inter
            print "All of %s and %s: %d" \
                %(fings[i], fings[j], len(inter))
    if len(fings) == 3:
        inter = bothStrinsUnion['F1'] & bothStrinsUnion['F2'] & bothStrinsUnion['F3']
        interUnion['F1', 'F2', 'F3'] = inter
        print "All of F1 F2 F3 : %d" %len(inter)

    # Unioning the union for both stringencies
    print 
    print "Unioning the unions"
    uUnion = {}
    for i in range(len(fings)):
        for j in range(i+1, len(fings)):
            union = bothStrinsUnion[fings[i]] | bothStrinsUnion[fings[j]]
            uUnion[fings[i], fings[j]] = union
            print "All of %s and %s: %d" \
                %(fings[i], fings[j], len(union))
    if len(fings) == 3:
        union = bothStrinsUnion['F1'] | bothStrinsUnion['F2'] | bothStrinsUnion['F3']
        uUnion['F1', 'F2', 'F3'] = union
        print "All of F1 F2 F3 : %d" %len(union)

    # Don't perform the next set unless we have the set of all proteins
    if allProts == None:
        return

    # For the canonical sets only.
    print
    print "How helpful are nearest neighbors for coverage?"
    print '#'*64
    for f in fings:
        for s in strins:
            missing = allProts - bindset[f,s]
            print "Missing from %s (%s)  : %d" %(f, s, len(missing))

            numHaveNeighbors1 = getNumHaveNeighbors(missing, bindset[f,s],
                                                    type = 'oneoff')
            print "Neighbor coverage %s %s (all one-off): %d" \
                %(f, s, numHaveNeighbors1 + len(bindset[f,s]))
            numHaveNeighbors2 = getNumHaveNeighbors(missing, bindset[f,s],
                                                    type = 'decomp-oneoff')
            print "Neighbor coverage %s %s (decomp one-off): %d" \
                %(f, s, numHaveNeighbors2 + len(bindset[f,s]))
            numHaveNeighbors3 = getNumHaveNeighbors(missing, bindset[f,s],
                                                    type = 'decomp-twooff-furthest')
            print "Neighbor coverage %s %s (decomp two-off furthest): %d" \
                %(f, s, numHaveNeighbors3 + len(bindset[f,s]))
            

    # How helpful for neighbors for intersected high and low stringencies?
    print
    print "UNIONING HIGH AND LOW STRINGENCIES"
    for k in sorted(bothStrinsUnion.keys()):
        missing = allProts - bothStrinsUnion[k]
        print "Missing from %s : %d" %(k, len(missing))
        numHaveNeighbors1 = getNumHaveNeighbors(missing, bothStrinsUnion[k],
                                                type = 'oneoff')
        print "Neighbor coverage %s (all one-off): %d" \
            %(k, numHaveNeighbors1 + len(bindset[f,s]))
        numHaveNeighbors2 = getNumHaveNeighbors(missing, bothStrinsUnion[k],
                                                type = 'decomp-oneoff')
        print "Neighbor coverage (decomp one-off): %d" \
            %(k, numHaveNeighbors2 + len(bindset[f,s]))
        numHaveNeighbors3 = getNumHaveNeighbors(missing, bothStrinsUnion[k],
                                                type = 'decomp-twooff-furthest')
        print "Neighbor coverage %s %s (decomp two-off furthest): %d" \
            %(k, numHaveNeighbors3 + len(bindset[f,s]))

    print "UNIONING STRINGENCIES AND FINGERS"
    for i in range(len(fings)):
        for j in range(i+1, len(fings)):
            twoFings = (bothStrinsUnion[i] | bothStrinsUnion[j])
            missing = allProts - twoFings
            print "Missing from %s %s: %d" %(fings[i], fings[j], len(missing))
            numHaveNeighbors = getNumHaveNeighbors(missing, twoFings,
                                                    type = 'oneoff')
            print "Neighbor coverage %s %s (all one-off): %d" \
                %(fings[i], fings[j], numHaveNeighbors + len(twoFings))
            numHaveNeighbors = getNumHaveNeighbors(missing, twoFings,
                                                    type = 'decomp-oneoff')
            print "Neighbor coverage %s %s (all one-off): %d" \
                %(fings[i], fings[j], numHaveNeighbors + len(twoFings))
            numHaveNeighbors = getNumHaveNeighbors(missing, twoFings,
                                                    type = 'decomp-twooff-furthest')
            print "Neighbor coverage %s %s (all one-off): %d" \
                %(fings[i], fings[j], numHaveNeighbors + len(twoFings))
    if len(fings) == 3:
        union = bothStrinsUnion['F1'] | bothStrinsUnion['F2'] | bothStrinsUnion['F3']
        missing = allProts - union
        print "Missing from F1 F2 F3: %d" %(len(missing))
        numHaveNeighbors = getNumHaveNeighbors(missing, union,
                                                type = 'oneoff')
        print "Neighbor coverage F1 F2 F3 (all one-off): %d" \
            %(numHaveNeighbors + len(union))
        numHaveNeighbors = getNumHaveNeighbors(missing, union,
                                                type = 'decomp-oneoff')
        print "Neighbor coverage F1 F2 F3 (all one-off): %d" \
            %(numHaveNeighbors + len(union))
        numHaveNeighbors = getNumHaveNeighbors(missing, union,
                                                type = 'decomp-twooff-furthest')
        print "Neighbor coverage F1 F2 F3 (all one-off): %d" \
            %(numHaveNeighbors + len(union))



    return bothStrinsUnion


def getNumHaveNeighbors(missing, bindset, type = 'oneoff'):
    # Given a set of missing canoncial sequences
    # and a set of binding sequences, returns the 
    # number of missing sequences that have 
    # a nearest neighbor in the set of binding seqs

    numHaveNeighbors = 0
    
    # Number have >= 1 neighbor only one off
    if type == 'oneoff':
        for seq in missing:
            neighbors = []
            for i in [0,1,2,3]:
                for a in aminos:
                    if a != seq[i]:
                        neighbors.append(seq[:i] + a\
                                         + seq[i+1:])
            for n in neighbors:
                if n in bindset:
                    numHaveNeighbors += 1
                    break

    elif type == 'decomp-oneoff':
        decomp = {1: [0,1,2], 2: [1,0,3], 3: [3,2,1]}
        
        for seq in missing:
            # Get neighbors for each base position
            for base in decomp.keys():
                neighbors = []
                for i in decomp[base]:
                    for a in aminos:
                        if a != seq[i]:
                            neighbors.append(seq[:i] + a\
                                             + seq[i+1:])

                # Check to see if any of this seqs neighbors are in 
                # the binding set
                if len(set(neighbors) & bindset) > 0:
                    neighborForEachBase = True
                else:
                    neighborForEachBase = False
                    break
                
            if neighborForEachBase:
                numHaveNeighbors += 1

    elif type == 'decomp-twooff-furthest':
        decomp = {1: [0,1,2], 2: [1,0,3], 3: [3,2,1]}

        for seq in missing:
            for base in decomp.keys():
                neighbors = []
                furthest = decomp[base][:2]
                closest = decomp[base][-1]

                # Vary two least important positions simultaneously
                for a1 in aminos:
                    for a2 in aminos:
                        newSeq = list(seq)
                        newSeq[furthest[0]] = a1
                        newSeq[furthest[1]] = a2
                        if ''.join(newSeq) != seq:
                            neighbors.append(''.join(newSeq))
                
                # Vary closest poition to fixed position by itself
                for a in aminos:
                    newSeq = list(seq)
                    newSeq[closest] = a
                    if ''.join(newSeq) != seq:
                            neighbors.append(''.join(newSeq))
                        
                # Check to see if any of this seqs neighbors are in 
                # the binding set
                if len(set(neighbors) & bindset) > 0:
                    neighborForEachBase = True
                else:
                    neighborForEachBase = False
                    break

            if neighborForEachBase:
                numHaveNeighbors += 1


    return numHaveNeighbors

def compareToNatural(bindset, natset, org, fings, 
                     strins, bothStrinsUnion):
	# Compares the binding set of to the natural set of 
	# ZFs.  bindset should be indexed by (fing,strin) tuples.
    print
    print "Fraction of %s ZFs (-1,2,3,6) Captured: (%d possible)" \
    	%(org, len(natset))
    print "#"*64
    for f in fings:
        for s in strins:
            inter = bindset[f,s] & natset
            print "%s %s: %.3f" %(f, s, len(inter)/float(len(natset)))
            haveNeighbors1 = getNumHaveNeighbors((natset - inter), bindset[f,s],    
                                                 'oneoff')
            haveNeighbors2 = getNumHaveNeighbors((natset - inter), bindset[f,s],    
                                                 'decomp-oneoff')
            haveNeighbors3 = getNumHaveNeighbors((natset - inter), bindset[f,s],    
                                                 'decomp-twooff-furthest')
            print "Neighbor Coverage (all one-off): %.3f" \
                %( (haveNeighbors1 + len(inter)) / float(len(natset)))
            print "Neighbor Coverage (decomp one-off): %.3f" \
                %( (haveNeighbors2 + len(inter)) / float(len(natset)))
            print "Neighbor Coverage (decomp two-off furthest): %.3f" \
                %( (haveNeighbors3 + len(inter)) / float(len(natset)))

    print
    for i in range(len(fings)):
        for j in range(i + 1, len(fings)):
            for s in strins:
                inter = (bindset[fings[i], s] & bindset[fings[j], s]) & natset
                union = (bindset[fings[i], s] | bindset[fings[j], s]) & natset
                print "%s %s Intersection (%s): %.3f" %(fings[i], fings[j], 
                                                        s, len(inter)/\
                                                        float(len(natset)))
                print "%s %s Union (%s)       : %.3f" %(fings[i], fings[j], 
                                                        s, len(union)/\
                                                        float(len(natset)))

    print            
    for s in strins:
    	inter = bindset[fings[0], s]
    	union = set()
    	for f in fings:
    		union = union | bindset[f,s]
    		inter = inter & bindset[f,s]
    	union = union & natset
        inter = inter & natset
        print "All Intersection (%s): %.3f" %(s, len(inter)/float(len(natset)))
        print "All Union (%s): %.3f" %(s, len(union)/float(len(natset)))

    print
    print "UNIONING HIGH AND LOW STRINGENCIES:"
    for k in sorted(bothStrinsUnion.keys()):
        inter = (bothStrinsUnion[k] & natset)
        print "All of %s: %.3f" %(k, len(inter)/float(len(natset)))

    for i, k1 in enumerate(sorted(bothStrinsUnion.keys())):
        for j, k2 in enumerate(sorted(bothStrinsUnion.keys())):
            if j > i:
                union = (bothStrinsUnion[k1] | bothStrinsUnion[k2]) & natset
                print "All of %s and %s: %.3f" \
                    %(k1, k2, len(union)/float(len(natset)))
                
                haveNeighbors1 = getNumHaveNeighbors((natset - union), bindset[f,s],    
                                                 'oneoff')
                haveNeighbors2 = getNumHaveNeighbors((natset - union), bindset[f,s],    
                                                     'decomp-oneoff')
                haveNeighbors3 = getNumHaveNeighbors((natset - union), bindset[f,s],    
                                                     'decomp-twooff-furthest')
                print "Neighbor Coverage (all one-off): %.3f" \
                    %( (haveNeighbors1 + len(union)) / float(len(natset)))
                print "Neighbor Coverage (decomp one-off): %.3f" \
                    %( (haveNeighbors2 + len(union)) / float(len(natset)))
                print "Neighbor Coverage (decomp two-off furthest): %.3f" \
                    %( (haveNeighbors3 + len(union)) / float(len(natset)))
                
    if len(fings) == 3:
        unionAll = (bothStrinsUnion['F1'] | bothStrinsUnion['F2'] \
                    | bothStrinsUnion['F3']) & natset
        print "All of F1, F2, and F3: %.3f" \
            %(len(unionAll)/float(len(natset)))

        haveNeighbors1 = getNumHaveNeighbors((natset - unionAll), bindset[f,s],    
                                         'oneoff')
        haveNeighbors2 = getNumHaveNeighbors((natset - unionAll), bindset[f,s],    
                                             'decomp-oneoff')
        haveNeighbors3 = getNumHaveNeighbors((natset - unionAll), bindset[f,s],    
                                             'decomp-twooff-furthest')
        print "Neighbor Coverage (all one-off): %.3f" \
            %( (haveNeighbors1 + len(unionAll)) / float(len(natset)))
        print "Neighbor Coverage (decomp one-off): %.3f" \
            %( (haveNeighbors2 + len(unionAll)) / float(len(natset)))
        print "Neighbor Coverage (decomp two-off furthest): %.3f" \
            %( (haveNeighbors3 + len(unionAll)) / float(len(natset)))


def compareBindingSets(bset1, bset2):
    inter = bset1 & bset2
    union = bset1 | bset2
    print "Set 1: %d" %len(bset1)
    print "Set 2: %d" %len(bset2)
    print "Intersection: %d" %len(inter)
    print "Union: %d" %len(union)
    print "Jaccard: %.3f" %(len(inter)/float(len(union)))

def computeBindingDiversity(proteinDir):
    # Get the protein sets from the binding data.

    sys.stdout = open('../stats/bindingDiv_' + \
                      '_'.join(proteinDir.split('_')[2:]) +\
                      '.txt', 'w')

    maxSize6 = 20**6
    maxSize4 = 20**4
    nucs = ['A', 'C', 'G', 'T']
    aminos = ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'H', 'K', 'L', 
              'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # Get the set of all possible canonical proteins
    allProts4 = []
    for i in product(aminos, aminos, aminos, aminos):
        allProts4.append(''.join(i))
    allProts4 = set(allProts4)

    # Get the dictionary of binding sets indexed by finger and stringency
    pathPref = '../data/b1hData/newDatabase/6varpos'
    #proteinDir = 'protein_seq_cut3bc_0_5'
    fings = ['F1', 'F2', 'F3']
    strins = ['high', 'low']
    bindingSets6 = {}
    bindingSets4 = {}
    for f in fings:
        for s in strins:
            bindingSets6[f,s] = getProtSet('/'.join([pathPref,f,s,proteinDir,'all.txt']), 
                                           range(6))
            bindingSets4[f,s] = getProtSet('/'.join([pathPref,f,s,proteinDir,'all.txt']),
                                           [0,2,3,5])

    printBindingSetStats(fings, strins, bindingSets6, maxSize6)
    bothStrinsUnion = printBindingSetStats(fings, strins, bindingSets4, 
                                           maxSize4, allProts4)

    dmelSet = getUniqueShilpaZFs('../data/shilpa/Drosophila_melanogaster_ZF.fulldom', 
                                 [0, 2, 3, 6])
    hsapSet = getUniqueShilpaZFs('../data/shilpa/Homo_sapiens_ZF.fulldom', 
                                 [0, 2, 3, 6])
    compareToNatural(bindingSets4, hsapSet, 'Human', fings, strins, bothStrinsUnion)
    compareToNatural(bindingSets4, dmelSet, 'Fly', fings, strins, bothStrinsUnion)

def main():

    protDirs = ['protein_seq_cut10bc_0', 'protein_seq_cut10bc_0_5', 
                'protein_seq_cut10bc_025', 'protein_seq_cut3bc_0', 
                'protein_seq_cut3bc_0_5', 'protein_seq_cut3bc_025']

    for protDir in protDirs:
        computeBindingDiversity(protDir)

if __name__ == '__main__':
	main()