#####################################
#
# Lorenz Chaotic Model 1963
# 
#####################################

#####################################
# Imports: 
import numpy as np
from scipy import integrate
#####################################

# Default: Lorenz chaotic model 1963
def lorenz_deriv((x, y, z), t0, sigma=10., beta=8./3, rho=28.0):
	# Compute the time-derivative of a Lorenz system.
	return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]

x0 = [0.1, 0, 0.1]  # starting vector
t = np.linspace(0, 1000, 100000)  # interval of .01


x_t = integrate.odeint(lorenz_deriv, x0, t)


print np.shape(x_t)


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(x_t[:,0], x_t[:,1], x_t[:,2])
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Lorenz Attractor")
plt.show()
