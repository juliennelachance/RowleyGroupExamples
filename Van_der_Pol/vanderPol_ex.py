#####################################
#
# van der Pol Oscillator
# 
#####################################

#####################################
# Imports: 
import numpy as np
from scipy import integrate
#####################################


def van_der_Pol_deriv((x, y), t0, mu=0.2):
	# Compute the time-derivative of a van der Pol system.
	return [y, -mu*(x**2.0 - 1.0)*y - x]

x0 = [-4.0, -4.0]  # starting vector
t = np.linspace(0, 100, 10000)  # interval of .01


x_t = integrate.odeint(van_der_Pol_deriv, x0, t)


print np.shape(x_t)


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca()

ax.plot(x_t[:,0], x_t[:,1])
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_title("van der Pol Oscillator")
plt.gca().set_aspect('equal')
plt.show()
