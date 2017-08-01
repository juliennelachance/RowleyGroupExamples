#####################################
#
# Von Karman Street Example
# From Marko Budisic's repository
#
#####################################



#####################################
# Imports: 
import numpy
from math import *
import pylab
import sys
import matplotlib.colors
import scipy.linalg
import scipy.io as sio
import numpy as np
#####################################

# -----------------------------------------------------
# Load the data and store as Data
# Data source: https://www.dropbox.com/s/p0cl8t7q9l2qwe4/VonKarmanStreet.mat?dl=1
# See also the GitHub of Marko Budisic: https://github.com/mbudisic/koopman/blob/f740d3accd6e3aaf822aeff5cdb4a7ad0ddae285/examples/VonKarmanStreet.m
Data = sio.loadmat('VonKarmanStreet.mat')  	# .mat file 
						# U contains a time evolution of velocity in format U( component, ycoord, xcoord, timestep )
						# Lattice-Boltzmann simulation of Von Karman street, simulated using flow by FlowKit. 
						# https://youtu.be/M2PqI2JD2jo
U = Data['U']					# Data matrix (4D) for analysis.
# -----------------------------------------------------

# Evaluate discrete derivatives used to compute vorticity
Ux = np.diff(U, 1, 2 );
Uy = np.diff(U, 1, 1 );

Uxadd = Ux[:,:,-1,:]
Uyadd = Uy[:,-1,:,:]

# Pad to preserve data size
Ux = np.append(Ux, Uxadd[:,:,np.newaxis,:], axis=2)
Uy = np.append(Uy, Uyadd[:,np.newaxis,:,:], axis=1)

# vorticity
Vorticity = Ux[1,:,:,:] - Uy[0,:,:,:];
print np.shape(Vorticity)

# speed
Speed = np.sqrt(sum(np.square(U),0));
print np.shape(Speed)

# stack data
Data = np.concatenate((Speed[np.newaxis,:,:,:], Vorticity[np.newaxis,:,:,:]), 0);
print np.shape(Data)
[nVars, nRows, nCols, nSteps] = np.shape(Data);


# Plotting:
step = 1

pylab.figure()
pylab.subplot(2,1,1)
pylab.pcolor(np.squeeze(Data[0,:,:,step]))

pylab.xlabel('X');
pylab.ylabel('Y');
pylab.title('Speed at step '+ repr(step))
pylab.axis('equal')
pylab.axis('tight')

pylab.subplot(2,1,2)
pylab.pcolor(np.squeeze(Data[1,:,:,step]))

pylab.xlabel('X');
pylab.ylabel('Y');
pylab.title('Vorticity at step '+ repr(step))
pylab.axis('equal')
pylab.axis('tight')
pylab.show()

























