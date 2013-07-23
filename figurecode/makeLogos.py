import os

def makeLogo(infile, outfile, format = 'pdf', alpha = 'protein',
             composition = 'equiprobable', size = 'small', 
             colScheme = 'chemistry', title = '', xlabel = '', 
             ylabel = ''):
	# Creates a logo based on a pwm by calling weblogo.
	#
	# - infile is a path to a pwm file in the transfac 
	# file format.
	# - outfile is the path where the file will be saved
	# The remaining parameters correspond ot flags for weblogo.
	# See weblogo CLI for more info about them.
	opts = '-F transfac'
	os.system('weblogo %s < %s > %s' %(opts, infile, outfile)) 
	pass