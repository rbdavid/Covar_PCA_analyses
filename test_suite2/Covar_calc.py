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
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

flatten = np.ndarray.flatten
zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
eigen = np.linalg.eig
flush = sys.stdout.flush

# ----------------------------------------
# VARIABLE DECLARATION
# ----------------------------------------

config_file = sys.argv[1]

necessary_parameters = ['pdb_file','traj_loc','start','end','average_pdb']
all_parameters = ['pdb_file','traj_loc','start','end','average_pdb','alignment','covar_selection','coarseness','fine_grain_selection','dist_covar_filename','cart_covar_filename','functionalize_bool','functionalized_dist_covar_filename','PCA_bool','PCA_eigenvalues_filename','PCA_eigenvectors_filename','summary_bool','summary_filename']

# ----------------------------------------
# SUBROUTINES:
# ----------------------------------------

def ffprint(string):
	print '%s' %(string)
	flush()

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''
	
	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['alignment'] = 'protein'
	parameters['covar_selection'] = 'protein'
	parameters['coarseness'] = 'COM'
	parameters['fine_grain_selection'] = None
	parameters['dist_covar_filename'] = 'distance_covar.dat'
	parameters['cart_covar_filename'] = 'cartesian_covar.dat'
	parameters['functionalize_bool'] = False
	parameters['functionalized_dist_covar_filename'] = 'functionalized_dist_covar.dat'
	parameters['PCA_bool'] = False
	parameters['PCA_eigenvalues_filename'] = 'PCA_eigenvalues_cart_covar.dat' 
	parameters['PCA_eigenvectors_filename'] = 'PCA_eigenvectors_cart_covar.dat' 
	parameters['summary_bool'] = True
	parameters['summary_filename'] = 'water_retention_analysis.summary'

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

	if parameters['coarseness'] not in ['COM','Atomic']:
		print "coarseness parameter does not match an acceptable value. Viable values are 'COM' and 'Atomic'. Killing job."
		sys.exit()

def summary(filename):
	with open(filename,'w') as W:
		W.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		W.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			W.write('%s ' %(sys.argv[i]))
		W.write('\n\nParameters used:\n')
		for i in all_parameters:
			W.write('%s = %s \n' %(i,parameters[i]))
		W.write('\n\n')

# ----------------------------------------
# MAIN:
# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# INITIATE THE AVG STRUCTURE; GRAB THE NECESSARY INFORMATION
ffprint('Initiating the average structure universe')
avg = MDAnalysis.Universe(parameters['average_pdb'])
avg_all = avg.select_atoms('all')
avg_align = avg.select_atoms(parameters['alignment'])
avg_all.translate(-avg_align.center_of_mass())
pos0 = avg_align.positions

# ----------------------------------------
# INITIALIZE THE ANALYSIS UNIVERSE; CREATE THE NECESSARY ATOM SELECTIONS
u = MDAnalysis.Universe(parameters['pdb_file'])
u_all = u.select_atoms('all')
u_align = u.select_atoms(parameters['alignment'])
u_covar = u.select_atoms(parameters['covar_selection'])

# ----------------------------------------
# COLLECTING NSTEPS DATA; ITERATE THROUGH ALL TRAJECTORIES AND SUM THE NUMBER OF TIMESTEPS
nSteps = 0
start = int(parameters['start'])
end = int(parameters['end'])
while start <= end:
#	u.load_new('%sproduction.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))
	u.load_new(parameters['traj_loc'])
	nSteps += len(u.trajectory)
	start += 1

ffprint('The number of timesteps to be analyzed is %d.' %(nSteps))

# ----------------------------------------
# COARSENESS -- Center Of Mass of Residues in the covar_selection

if parameters['coarseness'] == 'COM':
	ffprint('Performing a covariance analysis of the cartesian coordinates of the center of mass of residues defined in the covar_selection parameter.')

	nRes = u_covar.n_residues
	if nRes != avg.select_atoms(parameters['covar_selection']).n_residues:
		ffprint('The number of residues to be used in the covar_selection do not match between the average and analysis universes. Killing job.')
		sys.exit()

	# ----------------------------------------
	# MEMORY DECLARATION
	allCoord = zeros((nSteps,nRes,3),dtype=np.float64)
	avgCoord = zeros((nRes,3),dtype=np.float64)
	dist_covar_array = zeros((nRes,nRes),dtype=np.float64)
	msd_array = zeros(nRes,dtype=np.float64)

	# ----------------------------------------
	# TRAJECTORY ANALYSIS
	ffprint('Beginning trajectory analysis.')
	temp = 0 
	start = int(parameters['start'])
	while start <= end:
		ffprint('Loading trajectory %s' %(start))
#		u.load_new('%sproduction.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))
		u.load_new(parameters['traj_loc'])
		for ts in u.trajectory:
			u_all.translate(-u_align.center_of_mass())
			R,d = rotation_matrix(u_align.positions,pos0)		# MDAnalysis.analysis.align function
			u_all.rotate(R)
			for i in range(nRes):
				allCoord[temp,i,:] = u_covar.residues[i].center_of_mass()
				avgCoord[i,:] += u_covar.residues[i].center_of_mass()
			temp += 1
		start += 1
	avgCoord /= nSteps
	ffprint('Finished with the trajectory analysis. On to calculating the covariance matrix for residue-residue COM distance.')

	# ----------------------------------------
	# CALCULATING THE DISTANCE COVAR MATRIX OF RESIDUE-RESIDUE PAIRS 
	for i in range(nSteps):
		temp_array = zeros((nRes,3),dtype=np.float64)
		for res1 in range(nRes):
			temp_array[res1,:] = allCoord[i,res1,:] - avgCoord[res1,:]		# Calculating the delta r for every timestep
			msd_array[res1] += dot_prod(temp_array[res1,:],temp_array[res1,:])		# Sum over all timesteps; 
		for res1 in range(nRes):
			for res2 in range(res1,nRes):
				dist_covar_array[res1,res2] += dot_prod(temp_array[res1,:],temp_array[res2,:])
	dist_covar_array /= nSteps
	msd_array /= nSteps
	ffprint('Finished with filling the distance covariance matrix. On to normalizing the covar array.')
	
	# COMPLETE THE DISTANCE COVAR MATRIX ANALYSIS BY NORMALIZING BY THE VARIANCE
	for res1 in range(nRes):
		for res2 in range(res1,nRes):
			dist_covar_array[res1,res2] /= sqrt(msd_array[res1]*msd_array[res2])	# Normalizing each matrix element by the sqrt of the variance of the positions for res1 and res2
			dist_covar_array[res2,res1] = dist_covar_array[res1,res2]
	
	# ----------------------------------------
	# OUTPUTING THE DIST COVAR ARRAY
	with open(parameters['dist_covar_filename'],'w') as f:
		np.savetxt(f,dist_covar_array)

	# ----------------------------------------
	# CALCULATING THE CARTESIAN COVAR ARRAY OF RESIDUE RESIDUE PAIRS
	ffprint('Beginning the cartesian covariance matrix analysis')

	# ----------------------------------------
	# MEMORY DECLARATION
	cart_covar_array = zeros((3*nRes, 3*nRes),dtype=np.float64)
	cart_all_array = zeros(3*nRes,dtype=np.float64)
	cart_avg_array = zeros(3*nRes,dtype=np.float64)
	cart_msd_array = zeros(3*nRes,dtype=np.float64)

	# ----------------------------------------
	# CALCULATING THE CARTESIAN COVAR MATRIX OF RESIDUE-RESIDUE PAIRS 
	cart_avg_array = flatten(avgCoord)
	for i in range(nSteps):
		cart_all_array = flatten(allCoord[i])	# Each element in allCoord has three components (xyz); to get at the xyz components individually, flatten the first index
		temp_array = zeros(3*nRes,dtype=np.float64)
		for res1 in range(3*nRes):
			temp_array[res1] = cart_all_array[res1] - cart_avg_array[res1]
			cart_msd_array[res1] += temp_array[res1]*temp_array[res1]
		for res1 in range(3*nRes):
			for res2 in range(res1,3*nRes):
				cart_covar_array[res1,res2] += temp_array[res1]*temp_array[res2]
	
	cart_covar_array /= nSteps
	cart_msd_array /= nSteps
	
	# COMPLETE THE CARTESIAN COVAR MATRIX ANALYSIS BY NORMALIZING BY THE VARIANCE
	ffprint('Normalizing the cartesian covariance matrix using the variance.')
	for res1 in range(3*nRes):
		for res2 in range(res1,3*nRes):
			cart_covar_array[res1,res2] /= sqrt(cart_msd_array[res1]*cart_msd_array[res2])
			cart_covar_array[res2,res1] = cart_covar_array[res1,res2]	
	
	# OUTPUTING THE CARTESIAN COVAR ARRAY
	with open(parameters['cart_covar_filename'],'w') as f:
		np.savetxt(f,cart_covar_array)
	ffprint('Printed out the normalized cartesian covar array.')

# ----------------------------------------
# COARSENESS -- Atomic positions of fine_grain_selection in the covar_selection

elif parameters['coarseness'] == 'Atomic':
	ffprint('Performing a covariance analysis of the cartesian coordinates of the fine_grain_selection for the covar_selection.')

	u_fine_grain = u_covar.select_atoms(parameters['fine_grain_selection'])
	nAtoms = u_fine_grain.n_atoms

	avg_covar = avg.select_atoms(parameters['covar_selection'])

	if nAtoms != avg_covar.select_atoms(parameters['fine_grain_selection']).n_atoms:
		ffprint('The number of atoms to be analyzed in the fine_grain_selection of the covar_selection do not match between the average and analysis universes. Killing job.')
		sys.exit()

	# ----------------------------------------
	# MEMORY DECLARATION
	allCoord = zeros((nSteps,nAtoms,3),dtype=np.float64)
	avgCoord = zeros((nAtoms,3),dtype=np.float64)
	dist_covar_array = zeros((nAtoms,nAtoms),dtype=np.float64)
	msd_array = zeros(nAtoms,dtype=np.float64)

	# ----------------------------------------
	# TRAJECTORY ANALYSIS
	ffprint('Beginning trajectory analysis.')
	temp = 0 
	start = int(parameters['start'])
	while start <= end:
		ffprint('Loading trajectory %s' %(start))
#		u.load_new('%sproduction.%s/production.%s.dcd' %(parameters['traj_loc'],start,start))
		u.load_new(parameters['traj_loc'])
		for ts in u.trajectory:
			u_all.translate(-u_align.center_of_mass())
			R,d = rotation_matrix(u_align.positions,pos0)		# MDAnalysis.analysis.align function
			u_all.rotate(R)
			for i in range(nAtoms):
				allCoord[temp,i,:] = u_fine_grain.atoms[i].position
				avgCoord[i,:] += u_fine_grain.atoms[i].position
			temp += 1
		start += 1
	avgCoord /= nSteps
	ffprint('Finished with the trajectory analysis. On to calculating the covariance matrix for atom-atom positions.')

	# ----------------------------------------
	# CALCULATING THE DISTANCE COVAR MATRIX OF ATOM-ATOM PAIRS
	for i in range(nSteps):
		temp_array = zeros((nAtoms,3),dtype=np.float64)
		for atom1 in range(nAtoms):
			temp_array[atom1,:] = allCoord[i,atom1,:] - avgCoord[atom1,:]		# Calculating the delta r for every timestep
			msd_array[atom1] += dot_prod(temp_array[atom1,:],temp_array[atom1,:])		# Sum over all timesteps; 
		for atom1 in range(nAtoms):
			for atom2 in range(atom1,nAtoms):
				dist_covar_array[atom1,atom2] += dot_prod(temp_array[atom1,:],temp_array[atom2,:])
	dist_covar_array /= nSteps
	msd_array /= nSteps
	ffprint('Finished with filling the distance covariance matrix. On to normalizing the covar array.')
	
	# COMPLETE THE DISTANCE COVAR MATRIX ANALYSIS BY NORMALIZING BY THE VARIANCE
	for atom1 in range(nAtoms):
		for atom2 in range(atom1,nAtoms):
			dist_covar_array[atom1,atom2] /= sqrt(msd_array[atom1]*msd_array[atom2])	# Normalizing each matrix element by the sqrt of the variance of the positions for atom1 and atom2
			dist_covar_array[atom2,atom1] = dist_covar_array[atom1,atom2]
	
	# OUTPUTING THE DIST COVAR ARRAY
	with open(parameters['dist_covar_filename'],'w') as f:
		np.savetxt(f,dist_covar_array)

	# ----------------------------------------
	# CALCULATING THE CARTESIAN COVAR ARRAY OF ATOM-ATOM PAIRS
	ffprint('Beginning the cartesian covariance matrix analysis')

	# ----------------------------------------
	# MEMORY DECLARATION
	cart_covar_array = zeros((3*nAtoms, 3*nAtoms),dtype=np.float64)
	cart_all_array = zeros(3*nAtoms,dtype=np.float64)
	cart_avg_array = zeros(3*nAtoms,dtype=np.float64)
	cart_msd_array = zeros(3*nAtoms,dtype=np.float64)

	# ----------------------------------------
	# CALCULATING THE CARTESIAN COVAR MATRIX OF ATOM-ATOM PAIRS 
	cart_avg_array = flatten(avgCoord)
	for i in range(nSteps):
		cart_all_array = flatten(allCoord[i])	# Each element in allCoord has three components (xyz); to get at the xyz components individually, flatten the first index
		temp_array = zeros(3*nAtoms,dtype=np.float64)
		for atom1 in range(3*nAtoms):
			temp_array[atom1] = cart_all_array[atom1] - cart_avg_array[atom1]
			cart_msd_array[atom1] += temp_array[atom1]*temp_array[atom1]
		for atom1 in range(3*nAtoms):
			for atom2 in range(atom1,3*nAtoms):
				cart_covar_array[atom1,atom2] += temp_array[atom1]*temp_array[atom2]
	
	cart_covar_array /= nSteps
	cart_msd_array /= nSteps
	
	# COMPLETE THE CARTESIAN COVAR MATRIX ANALYSIS BY NORMALIZING BY THE VARIANCE
	ffprint('Normalizing the cartesian covariance matrix using the variance.')
	for atom1 in range(3*nAtoms):
		for atom2 in range(atom1,3*nAtoms):
			cart_covar_array[atom1,atom2] /= sqrt(cart_msd_array[atom1]*cart_msd_array[atom2])
			cart_covar_array[atom2,atom1] = cart_covar_array[atom1,atom2]	
	
	# OUTPUTING THE CARTESIAN COVAR ARRAY
	with open(parameters['cart_covar_filename'],'w') as f:
		np.savetxt(f,cart_covar_array)

# ----------------------------------------
# FUNCTIONALIZING THE DISTANCE CORRELATION MATRIX (FOR USE IN WISP AND VISUALIZATION)
if parameters['functionalize_bool']:
	ffprint('Beginning to functionalize the distance covar matrix.')
	nElements = len(dist_covar_array)
	for element1 in range(nElements):
		for element2 in range(element1,nElements):
			dist_covar_array[element1,element2] = -np.log(np.fabs(dist_covar_array[element1,element2]))
			dist_covar_array[element2,element1] = dist_covar_array[element1,element2]
	
	# OUTPUTING THE DIST COVAR ARRAY
	with open(parameters['functionalized_dist_covar_filename'],'w') as f:
		np.savetxt(f,dist_covar_array)
else: 
	ffprint('Functionalize != True. Not functionalizing (taking -log(|C_ij|)) the distance covar matrix.')

# ----------------------------------------
# PCA ANALYSIS OF CARTESIAN COVAR ARRAY
if parameters['PCA_bool']:
	ffprint('Beginning PCA analysis of the Cartesian covariance matrix.')
	eigval,eigvec = eigen(cart_covar_array)
	idx = eigval.argsort()[::-1]
	eigval = eigval[idx]
	
	nVec = len(eigvec)
	cumulative_eigval = zeros(nVec,dtype=np.float64)
	total_eigval = 0
	for i in range(nVec):
		total_eigval += eigval[i]
		cumulative_eigval[i] = total_eigval
	
	with open(parameters['PCA_eigenvalues_filename'],'w') as f:
		for i in range(nVec):
			f.write('%f   %f   %f   %f\n' %(eigval[i],eigval[i]/total_eigval,cumulative_eigval[i],cumulative_eigval[i]/total_eigval))
	
	with open(parameters['PCA_eigenvectors_filename'],'w') as f:
		for i in range(nVec):
			for j in range(nVec):
				f.write('%f   ' %(eigvec[j,i]))		# Writing each vector on one row/line now, instead of the vectors corresponding to columns in the eigvec array...; NOT projecting covar array onto the eigenvectors (do so outside of this damn script)
			f.write('\n')
else:
	ffprint('PCA_bool != True. Not performing a PCA analysis on the Cartesian Covar matrix.')

