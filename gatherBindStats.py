# Basic analysis of the diverity of the binding proteins.

import os
import sys
import math

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

def getProtDict(path, contacts, cut = 0):
	# Returns a dictionary mapping unique proteins 
	# (when using only given list of contact positions)
	# to a list of tuples (3-mer target bound, count).
	# path should be to a file the 'all.txt' format
	# and contatcs positions in file are -1,1,2,3,5,6
	# Can optionally remove proteins for which all
	# binding counts are less than the cut parameter.
	
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
		if freq >= cut:
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

def printBindingSetStats(fings, strins, bindset, maxSize):
    numpos = math.log(maxSize, 20)
    print "Max Size for %d positions: %d" %(numpos, maxSize)
    print "#"*64
    for f in fings:
        for s in strins:
            setSize = len(bindset[f,s])
            print "Size of %s %s: %d" %(f, s, setSize) 
        
    print
    print "Is high a subset of low? (|High - Low|/|High|)"
    print "#"*64
    for f in fings:
        highMinLow = bindset[f,'high'] - bindset[f,'low']
        print "%s: %.3f" %(f, len(highMinLow)/ \
                           float(len(bindset[f,'high'])))
    
    print 
    print "Combining Sets"
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


def compareToNatural(bindset, natset, org, fings, strins):
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