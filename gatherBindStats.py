# Basic analysis of the diverity of the binding proteins.

import os
import sys
import math
import re

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
    # Of those sequences not contained in a given set, how many 
    # have a nearest neighbor in the set if we vary positions
    # (-1, 3, 6)?
    print
    print "How helpful are nearest neighbors for coverage? (-1, 3, 6)"
    print '#'*64
    for f in fings:
        for s in strins:
            missing = allProts - bindset[f,s]
            numHaveNeighbors = getNumHaveNeighbors(missing, bindset[f,s])
            print "Missing from %s (%s)  : %d" %(f, s, len(missing))
            print "Have neighbors %s (%s): %d" \
                %(f, s, numHaveNeighbors)
            print "Total coverage %s %s (including neighbors): %d" \
                %(f, s, numHaveNeighbors + len(bindset[f,s]))

    # How helpful for neighbors for intersected high and low stringencies?
    print
    print "UNIONING HIGH AND LOW STRINGENCIES"
    for k in sorted(bothStrinsUnion.keys()):
        missing = allProts - bothStrinsUnion[k]
        numHaveNeighbors = getNumHaveNeighbors(missing, bothStrinsUnion[k])
        print "Missing from %s: %d"   %(k, len(missing))
        print "Have neighbors %s: %d" %(k, numHaveNeighbors)
        print "Total coverage %s (including neighbors): %d" \
            %(k, numHaveNeighbors + len(bothStrinsUnion[k]))

    return bothStrinsUnion


def getNumHaveNeighbors(missing, bindset):
    # Given a set of missing canoncial sequences
    # and a set of binding sequences, returns the 
    # number of missing sequences that have 
    # a nearest neighbor in the set of binding seqs

    numHaveNeighbors = 0
    for seq in missing:
        neighbors = []
        for i in [0,2,3]:
            for a in aminos:
                if a != seq[i]:
                    neighbors.append(seq[:i] + a\
                                     + seq[i+1:])
        for n in neighbors:
            if n in bindset:
                numHaveNeighbors += 1
                break

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
    if len(fings) == 3:
        unionAll = (bothStrinsUnion['F1'] | bothStrinsUnion['F2'] \
                    | bothStrinsUnion['F3']) & natset
        print "All of F1, F2, and F3: %.3f" \
            %(len(unionAll)/float(len(natset)))


def compareBindingSets(bset1, bset2):
    inter = bset1 & bset2
    union = bset1 | bset2
    print "Set 1: %d" %len(bset1)
    print "Set 2: %d" %len(bset2)
    print "Intersection: %d" %len(inter)
    print "Union: %d" %len(union)
    print "Jaccard: %.3f" %(len(inter)/float(len(union)))

def main():

	path = "../data/shilpa/Drosophila_melanogaster_ZF.fulldom"
	flyZFs = getUniqueShilpaZFs(path, [0, 1, 2, 5])
	print len(flyZFs)

if __name__ == '__main__':
	main()