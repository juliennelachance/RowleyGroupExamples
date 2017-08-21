#####################################
#
# Chua's Circuit
#
# Formulation courtesy of https://github.com/brendabrandy/chuacircuit/blob/master/chuacircuit/diffsol.py
# 
#####################################

#####################################
# Imports: 
import numpy as np
from scipy import integrate
from cmath import sin, cos, exp, pi, log, polar, rect, phase, sqrt
#####################################


# Parameters: 
c1 = 1.0 / 9.0
c2 = 1.0 
G = 0.7
L = 1.0 / 7.0 


# Initial conditions: 
x0 = [0.1, 0.1, 0.0001]  # starting vector

# Chua Circuit: 3-segment piecewise-linear function
def f(v):
	if v > -1 and v < 1:	# if v is between -1 and 1
		return -0.8* v
	elif v > 1: 		# if v is bigger than 1
		return -0.8 - 0.5*(v-1.0)
	else: 			# if v is less than -1
		return (0.8 - 0.5*(v+1.0))
# Chua Circuit: derivative function:
def chua_deriv((x, y, z), t):
	dx = np.empty([3])
	# Compute the time-derivative.
	dx[0] = (G*(y-x) - f(x) ) / c1	#dx/dt
	dx[1] = (G*(x - y) + z) / c2	#dy/dt
	dx[2] = - y / L			#dz/dt
	return dx



# Creating the time array:
t_start = 0.0				# initial time
t_end = 300.0				# final time
N = 100000				# number of steps
h = (t_end-t_start) / float(N)        	# step size
t = np.arange(t_start,t_end,h) 	# time array   

# Performing integration:
x_t = integrate.odeint(chua_deriv, x0, t)

# Plotting: 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
#ax = fig.gca()

ax.plot(x_t[:,0], x_t[:,1], x_t[:,2])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Chua Circuit: Phase Diagram")
plt.show()

