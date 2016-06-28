#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
##!/mnt/lustre_fs/users/mjmcc/apps/python2.7/bin/python
# ----------------------------------------
# USAGE:


# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

pdb_file = sys.argv[1]
traj_loc = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
avg_pdb = sys.argv[5]

alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'
important = 'protein'

zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

# ----------------------------------------
# MAIN PROGRAM:

# INITIATE AND CALCULATE THE IMPORTANT INFORMATION ABOUT THE AVG STRUCTURE
avg = MDAnalysis.Universe(avg_pdb)
avg_all = avg.select_atoms('all')
avg_align = avg.select_atoms(alignment)
avg_important = avg.select_atoms(important)

avg_all.translate(-avg_align.center_of_mass())
pos0 = avg_align.positions

nRes = len(avg_important.residues)
avgCoord = zeros((nRes,3),dtype=float)

for i in range(nRes):
	res0 = avg_important.residues[i]
	com0 = res0.center_of_mass()
	avgCoord[i,:] = com0

# INITIATE AND CREATE THE IMPORTANT ATOM SELECTIONS FOR THE IMPORTANT UNIVERSE
u = MDAnalaysis.Universe(pdb_file)
u_all = u.select_atoms('all')
u_align = u.select_atoms(alignment)
u_important = u.select_atoms(important)

if nRes != len(u_important.residues):
	ffprint('Number of residues do not match between average structure and trajectory universes')
	sys.exit()

# COLLECTING NSTEPS DATA... NEED TO ITERATE THROUGH ALL TRAJECTORIES TO COLLECT THIS...
nSteps = 0
temp = start
while temp <= end:
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,temp,temp))
	nSteps += len(u.trajectory)
	temp += 1

# ARRAY DECLARATION
allCoord = zeros((nSteps,nRes,3),dtype=float)
covar_array = zeros((3*nRes, 3*nRes),dtype=float)
temp_array = zeros(3*nRes,dtype=float)

# TRAJECTORY ANALYSIS; COLLECTING THE COM INFO OF EACH RESIDUE OF INTEREST
temp = 0
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))

	for ts in u.trajectory:
		dimensions = u.dimensions[:3]
		u_all.translate(-u_align.center_of_mass())

		R,d = rotation_matrix(u_align.positions,pos0)
		u_all.rotate(R)

		for i in range(nRes):
			res0 = u_important.residues[i]
			com0 = res0.center_of_mass()
			allCoord[temp,i,:] = com0

		temp += 1 
	start += 1

# CALCULATING THE COVAR ARRAY OF COM CARTESIAN COORDINATES
for i in range(nSteps):
	temp_array = flatten(allCoord[i])
	for res1 in range(3*nRes):
		for res2 in range(res1,3*nRes):
			covar_array[res1,res2] += temp_array[res1]*temp_array[res2]

covar_array /= nSteps

# COMPLETE THE COVAR MATRIX ANALYSIS; NEED TO SUBTRACT OUT THE MEAN
temp_array = flatten(avgCoord)
for res1 in range(3*nRes):
	for res2 in range(res1,3*nRes):
		covar_array[res1,res2] -= temp_array[res1]*temp_array[res2]
		covar_array[res2,res1] = covar_array[res1,res2]

# OUTPUTING THE COVAR MATRIX
out1 = open('%03d.%03d.cover_matrix.dat' %(sys.argv[3],end),'w')
for i in range(3*nRes):
	for j in range(3*nRes):
		out1.write('%e   ' %(covar_array[i,j]))
	out1.write('\n')
out1.close()

# WRITE A PCA ANALYSIS OF THIS COVAR ARRAY? 








