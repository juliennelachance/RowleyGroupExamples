#####################################
#
# POD with sklearn 
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
import modred as mr
import matplotlib.pyplot as plt
from cmath import sin, cos, exp, pi, log, polar, rect, phase, sqrt
from future.builtins import range
import scipy.ndimage.filters
from scipy import signal
import copy
from sklearn.decomposition import PCA
import numpy as np
from numpy.testing import assert_array_almost_equal
import time 
#####################################
t = time.time()



# User-defined option: 
n_modes = 2  # Should be >40 for full analysis... check eigenvalue semilog plot



# -----------------------------------------------------
# Dependencies: 
def do_sklearn_pca(X):
	X = X.T
	pca = PCA(n_components=n_modes)
	pca.fit(X)

	X_ = pca.transform(X)
	X_projected = pca.inverse_transform(X_)
	X_ = X_.T
	X_projected = X_projected.T

	loss = ((X_true[:,:len(X_projected[0,:])] - X_projected[:len(X_true),:]) ** 2).mean()
	print loss
	return X_, X_projected, pca

def do_eigenvalue_plots(X_in):
	U1, s1, V1h = numpy.linalg.svd(X_in, full_matrices=False)
	# Plot to determine reasonable truncation mode values: 
	pylab.semilogy(s1, '-x')
	pylab.title('Sorted eigenvalue plot of X')
	pylab.xlabel('Element index')
	pylab.ylabel('eigenvalue')
	pylab.show()

def show_pca_plots_clip(X_projected, start_val, end_val):
	N_c = 20
	v_min = -10
	v_max = 500
	min_col = min(len(X_projected[0,start_val:end_val]), len(X_true[0,start_val:end_val]))
	pylab.figure()
	pylab.subplot(3, 1, 1)
	pylab.contourf(np.clip(X_projected[:50,:(end_val-start_val)], v_min, v_max), N_c, vmin = v_min, vmax = v_max)
	pylab.colorbar()
	pylab.xticks( np.linspace(0, (end_val-start_val), 10),   np.around(np.linspace(start_val+1, end_val+1, 10),0))
	#pylab.xlabel('Carrington Rotation Number')
	pylab.ylabel('Sine(Latitude)')
	pylab.title('X Projected')
	pylab.subplot(3, 1, 2)
	pylab.contourf(np.clip(X_true[:50,start_val:end_val], v_min, v_max), N_c, vmin = v_min, vmax = v_max)
	pylab.xticks( np.linspace(0, (end_val-start_val), 10),   np.around(np.linspace(start_val+1, end_val+1, 10),0))
	#pylab.xlabel('Carrington Rotation Number')
	pylab.ylabel('Sine(Latitude)')
	pylab.colorbar()
	pylab.title('X Original')

	v_min_diff = -50
	v_max_diff = 50

	pylab.subplot(3, 1, 3)
	pylab.contourf( np.clip( (X_true[:50,start_val:end_val]) - (X_projected[:50,:(end_val-start_val)]), v_min_diff, v_max_diff), N_c)
	pylab.xticks( np.linspace(0, (end_val-start_val), 10),   np.around(np.linspace(start_val+1, end_val+1, 10),0))
	pylab.xlabel('Carrington Rotation Number')
	pylab.ylabel('Sine(Latitude)')
	pylab.colorbar()
	pylab.title('Difference')
	pylab.show()

def plot_reconstructed_spatial_means(plot_ans):
	pylab.figure()

	avg_arr_orig = np.mean(X_true[:,start_val:end_val], axis=0)
	pylab.plot(avg_arr_orig)

	#avg_arr_re = np.mean(plot_ans[:50,start_val:end_val], axis=0)
	avg_arr_re = np.mean(plot_ans, axis=0)
	limval = len(plot_ans[0,:])+1
	pylab.plot(avg_arr_re, 'r')
	pylab.xticks( np.linspace(0, (end_val-start_val), 10),   np.around(np.linspace(start_val, end_val, 10),0))
	pylab.title('Reconstructed Spatial Means')
	pylab.xlabel('Carrington Rotation Number')
	pylab.ylabel('Mean Sunspot Number Per Carr. Rot. Num.')
	pylab.legend(['Exact', 'DMD'])

	pylab.xlim(0,limval)
	pylab.show()
	return

def make_hankel(X_in):
	overlaps=160;
	print('The overlaps value is: ')
	print overlaps 
	H=numpy.zeros((len(X_in)*overlaps,len(X_in[0,:])-overlaps))
	for i in range(0,len(X_in[0,:])-overlaps):
		for j in range(0,overlaps):
			H[len(X_in)*(j):(len(X_in)*(j+1)),i]=X_in[:,(i+j)]
	return H

def run_dmd(X_in):
	X1 = X_in[:,:-1]
	Y1 = X_in[:,1:]
	U,s,V = numpy.linalg.svd(X1, full_matrices=False)
	modes = len(V) #can be truncated
	s = s[:modes]
	S = numpy.zeros((modes, modes), dtype=complex)
	S[:len(s), :len(s)] = numpy.diag(1.0/s)
	Xpinv=numpy.dot(V.T, numpy.dot(S,U.T))
	A = numpy.dot(Y1, Xpinv)
	v,d = sorted_eig(A)
	return A, v, d

def scale_vector(scale, vector):
	result = [0]*len(vector)
	for i in range(len(result)):
		result[i] = scale * vector[i]
	return result
def real_vector(vector):
	return map(lambda x: -1.0*x.real, vector)
def imag_vector(vector):
	return map(lambda x: x.imag, vector)
def sorted_eig(M, left=False):
	# Function to sort eigenvalues and normalize for POD
	w,v = scipy.linalg.eig(M, left=left, right=(not left))
	v = v.T
	sorted_wv = sorted([(abs(w[i]), list(v[i]/numpy.linalg.norm(v[i])), w[i]) for i in range(len(w))], reverse=True)
	_,v,w = zip(*sorted_wv)
	return numpy.array(w), numpy.array(v).T
def on_pick(event):
        ind = event.ind
	index = numpy.take(ind, 0)
	real_eig = numpy.take(numpy.real(numpy.take(v, ind)),0)
	imag_eig = numpy.take(numpy.imag(numpy.take(v, ind)),0)
        print('Eigenvalue Index: ' + repr(index))
	print('Value: ' + repr(numpy.take(numpy.take(v, ind),0)))
	imaglogeigv = numpy.imag(numpy.log(numpy.take(numpy.take(v, ind),0)))
	per = float(2.0*pi/imaglogeigv*27.2753/365)
        print('Approx. Period in years: ' + repr(per))
	print
	#do_mode_plot(index)

def plot_complex_eigenvalues():
	# Plot complex eigenvalues: 
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(numpy.real(v), numpy.imag(v), picker=True)
	ax.set_title('Eigenvalues of A')
	ax.set_xlabel('Real part of eigenvalues')
	ax.set_ylabel('Imag part of eigenvalues')
	# Generate numbers around the complex unit circle.
	N = 128
	theta = scale_vector(2*pi/N, range(N))
	exp_theta = map(lambda x: exp(1j * x), theta)
	real_part = real_vector(exp_theta)
	imag_part = imag_vector(exp_theta)
	pylab.plot(real_part,imag_part)
	fig.canvas.callbacks.connect('pick_event', on_pick)
	pylab.show()

def propagate_solution(start_val, end_val, X_prop, index):
	ans_=numpy.empty((len(X_prop),(end_val-start_val)))
	tmp_v = X_prop[:,index]
	ans_[:,0] = np.real(tmp_v)
	for i in range(1,(end_val-start_val)):
		ans_[:,i] = np.real(tmp_v)
		#tmp_v = ans_[:,i]
		tmp_v = numpy.dot(A, tmp_v)
	return ans_
# -----------------------------------------------------


# -----------------------------------------------------
# Load the data and store in X: 
# Data source: https://solarscience.msfc.nasa.gov/greenwch/bflydata.txt
mat_contents = sio.loadmat('sunspot_data.mat')  # MatLab file of butterfly diagram data
X = mat_contents['data_v']			# Data matrix (50x1907) for analysis.
X_true = copy.copy(X)				# A copy of the exact data for comparison later. 
# -----------------------------------------------------

# Determine how many PCA modes to keep: 
if(0==1):
	do_eigenvalue_plots(X)

# Perform PCA projection:
X_, X_projected, pca = do_sklearn_pca(X)

# Visualize how well PCA reconstruction fits the original data: 
start_val = 1
end_val = 1907
if(1==1):
	show_pca_plots_clip(X_projected, start_val, end_val)
if(1==1):
	plot_reconstructed_spatial_means(X_projected)

# Make Hankel matrix:
X_ = make_hankel(X_)

# DMD application (with 160 time delays):
A, v, d = run_dmd(X_)

# Visualize eigenvalues. Click plot to see corresponding period.
if(1==1):
	plot_complex_eigenvalues()

# Further analysis: Can reconstruct using A matrix from DMD.


