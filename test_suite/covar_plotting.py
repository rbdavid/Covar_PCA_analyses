#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

from plotting_functions import *
import sys

dist_data_file = sys.argv[1]
func_data_file = sys.argv[2]

dist_data = np.loadtxt(dist_data_file)
func_data = np.loadtxt(func_data_file)

matrix2d(dist_data,'Residue Number','Residue Number', 'Covariance','Covar','test_system') #,vmin=-1.00,vmax=1.00)
matrix2d(func_data,'Residue Number','Residue Number', 'Covariance','func_Covar','test_system') #,vmin=-1.00,vmax=1.00)

