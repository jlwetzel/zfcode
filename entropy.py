import numpy as np

def entropy(f):
	# Returns the shannon entropy of the distribution 
	# f using base 2 logs and using the convention that 
	# plg(p) = 0 when p = 0
	#
	# f should be array-like

	ent = 0.0
	for i in f:
		if i == 0:
			continue
		else:
			ent -= i * np.log2(i)

	return ent