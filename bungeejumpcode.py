"""
Program: Dynamical Systems Template
File: dynamicalSystemsTemplate.pylab
Author: blakekas
Version: 2-10-2018.1
Summary:
    Basic script for solving a system of first order differential equations
    that describes a dynamical system.
Usage: python dynamicalSystemsTemplate.py
Version History:
    1-15-2018.1: base
    12-11-2018.1: added some documentation
    
"""

import scipy.integrate as si
import numpy as np
import matplotlib.pylab as plt

def step(x):
    return 1 * (x > 0)

# Create a function that defines the rhs of the differential equation system
def calc_RHS(y,t,p):
#   y = a list that contains the system state
#   t = the time for which the right-hand-side of the system equations
#       is to be calculated.
#   p = a tuple that contains any parameters needed for the model
#
    import numpy as np

#   Unpack the state of the system
    y0 = y[0]
    y1 = y[1]

#   Unpack the parameter list
    m , L , k1 , k2 , X1 , kd , = p

#   Calculate the bungee force
    if y0 < L:
        Fb = 0
    elif y0-L < X1:  
        Fb = -k1*(y0-L)
    else:
        Fb= -k2*(y0-L)

#   Calculate the rates of change (the derivatives)
    dy0dt = y1
    dy1dt = 9.81 + Fb/m -kd*y1**2*np.sign(y1)

    return [dy0dt,dy1dt]

# Define the initial conditions
y_0 = [0.,0.]

# Define the time grid
t = np.linspace(0,120,200)

# Define the model parameters
m = 75
L = 3.
k2 = 162.
k1 = 77.
X1 = 4.88
kd = 1.e-2
p = m , L , k1 , k2 , X1 , kd ,

# Solve the DE
sol = si.odeint(calc_RHS,y_0,t,args=(p,))
y0 = sol[:,0]
y1 = sol[:,1]

# Plot the solution
plt.plot(t,y0,color='g')
#plt.plot(t,y1,color='b')
plt.show()
