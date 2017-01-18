#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
# ----------------------------------------
# USAGE:
# ----------------------------------------


# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import numpy as np
try:
	import cPickle as pickle
except:
	import pickle

# ----------------------------------------
# VARIABLE DECLARATION
# ----------------------------------------

data_file = sys.argv[1]		# Assumes data file is a ascii file of floats (rows by columns)
output_file = sys.argv[2]

flush = sys.stdout.flush
test = True

# ----------------------------------------
# SUBROUTINES:
# ----------------------------------------

def ffprint(string):
	print '%s' %(string)
	flush()

# ----------------------------------------
# MAIN PROGRAM:
# ----------------------------------------

data1 = np.loadtxt(data_file)
pickle.dump(data1,open('%s' %(output_file),'wb'))

if test == True:
	data2 = pickle.load(open('%s' %(output_file),'rb'))
	print 'Equal?:', (data1 == data2)

