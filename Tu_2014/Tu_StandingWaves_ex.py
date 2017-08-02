#####################################
#
# Standing Waves Example
# Tu et al., 2014 
# 
#####################################

#####################################
# Imports: 
import numpy as np
from math import *
import matplotlib.pyplot as plt
#####################################

def van_der_Pol_deriv((x, y), t0, mu=0.2):
	# Compute the time-derivative of a van der Pol system.
	return [y, -mu*(x**2.0 - 1.0)*y - x]

q = 1
uv0 = [q, 0]  # starting vector
theta = pi/6
t = np.linspace(0, 0.001, 50)  # interval of .01

u = np.zeros(len(t))
v = np.zeros(len(t))
u[0] = uv0[0]
v[0] = uv0[1]
for i in range(1, len(t)):
	u[i] = (cos(theta))*u[i-1] - (sin(theta))*v[i-1]
	v[i] = (sin(theta))*u[i-1] + (cos(theta))*v[i-1]

fig = plt.figure()
ax = fig.gca()
plt.plot(u, 'b', label='u')
plt.plot(v, 'r', label='v')
plt.legend()
ax.set_xlabel("Time")
ax.set_ylabel("Amplitude")
ax.set_title("Standing Waves")
plt.show()
