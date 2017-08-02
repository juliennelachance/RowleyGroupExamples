#####################################
#
# Stochastic Dynamics Example
# Tu et al., 2014 
# 
#####################################

#####################################
# Imports: 
import numpy as np
from math import *
import matplotlib.pyplot as plt
#####################################

# z[k+1] = lambda*z[k] + n[k], where
# 	lambda = 0.5 and n[k] is white noise with variance sig^2 = 10

num_samples = 1000
lambda_v = 0.5
nk = np.random.normal(0, sqrt(10.0), size=num_samples)

zk = np.zeros(num_samples)
for i in range(1,num_samples):
	zk[i] = lambda_v*zk[i-1] + nk[i-1]

fig = plt.figure()
ax = fig.gca()
plt.plot(zk)
ax.set_xlabel("k")
ax.set_ylabel("zk")
ax.set_title("Typical Trajectory for Zero I.C.")
plt.show()

fig = plt.figure()
ax = fig.gca()
plt.scatter(zk[0:-2], zk[1:-1])
plt.gca().set_aspect('equal')
ax.set_xlabel("z[k]")
ax.set_ylabel("z[k+1]")
ax.set_title("Correlation of z[k] and z[k+1]")
plt.plot([-20, 20], [-10, 10], 'r-')
plt.show()
