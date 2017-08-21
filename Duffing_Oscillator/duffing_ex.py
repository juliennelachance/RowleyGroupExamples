#####################################
#
# Duffing Oscillator
# 
#####################################

#####################################
# Imports: 
import numpy as np
from scipy import integrate
from cmath import sin, cos, exp, pi, log, polar, rect, phase, sqrt
#####################################

unforced_ex = 1		# Set as 0 for example without forcing; 1 for example with forcing. 

# Parameters: 
omega = 1.0
gamma = 0.1
epsilon = 0.25
Omega = 2.0
if(unforced_ex == 0):
	Gamma = 0.0
else:
	Gamma = 2.5

# Initial conditions: 
if(unforced_ex == 0):
	x0 = [3.0, 10.0]  # starting vector
else:
	x0 = [0.0, 0.0]



# Duffing Oscillator: derivative function:
def duff_deriv((x, y), t):
	dx = np.empty([2])
	# Compute the time-derivative.
	dx[0] = y;
	dx[1] = -gamma*y+omega*x - epsilon*np.power(x,3) + Gamma*cos(Omega*t);
	return dx

# Creating the time array:
if(1==0):
	h = 1e-1 	# time step
	period = 2*np.pi/(1.0*Omega)/100.0
	T = 40000 	# length of the simulation
	t = np.arange(0,T,period)
else:
	t = np.arange(0,100,(100.0/5000.0))

# Performing integration:
x_t = integrate.odeint(duff_deriv, x0, t)

# Plotting: 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
#ax = fig.gca(projection='3d')
ax = fig.gca()

ax.plot(x_t[:,0], x_t[:,1])
ax.set_xlabel("x")
ax.set_ylabel("x dot")
ax.set_title("Duffing Oscillator: Phase Diagram")
plt.show()








