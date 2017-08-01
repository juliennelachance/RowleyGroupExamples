########################################################
#
# The following code is courtesy of Hao Zhang.
#
# "An example to demonstrate online dynamic mode decomposition
# 
# We take a 2D time varying system given by dx/dt = A(t)x
# where x = [x1,x2]', A(t) = [0,w(t);-w(t),0], 
# w(t)=1+epsilon*t, epsilon=0.1. The slowly time varying eigenvlaues of A(t)
# are pure imaginary, +(1+0.1t)j and -(1+0.1t)j, where j is the imaginary unit."
#
########################################################

import numpy as np
from scipy.integrate import odeint
import time
import matplotlib.pyplot as plt


# define dynamics
epsilon = 1e-1
def dyn(x,t):
    x1, x2 = x
    dxdt = [(1+epsilon*t)*x2,-(1+epsilon*t)*x1]
    return dxdt
# integrate from initial condition [1,0]    
tspan = np.linspace(0,10,101)
dt = 0.1
x0 = [1,0]
xsol = odeint(dyn,x0,tspan).T
# extract snapshots
x, y = xsol[:,:-1], xsol[:,1:]
t = tspan[1:]
# true dynamics, true eigenvalues
n, m = len(x[:,0]), len(x[0,:])
A = np.empty((n,n,m))
evals = np.empty((n,m),dtype=complex)
for k in range(m):
    A[:,:,k] = np.array([[0,(1+epsilon*t[k])],[-(1+epsilon*t[k]),0]])
    evals[:,k] = np.linalg.eigvals(A[:,:,k])


# visualize snapshots
plt.figure()
plt.plot(tspan, xsol[0,:], 'bs-', linewidth=2.0,  label='$x_1(t)$')
plt.plot(tspan, xsol[1,:], 'g^-', linewidth=2.0,  label='$x_2(t)$')
plt.legend(loc='best',fontsize=20 ,shadow=True)
plt.xlabel('Time', fontsize=20)
plt.title('Snapshots', fontsize=20)
plt.tick_params(labelsize=20)
plt.grid()
plt.show()
