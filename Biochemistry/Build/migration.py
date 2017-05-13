# biochemical signaling model 

# import pydelay and numpy and pylab

from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from pydelay import dde23
import csv
import json

k1_val = 0.001
qR_val = 0.1

# the model equations 
eqns = { 
    'Ac' : '-k1*R*Ac',
    'A'  : 'k1*R(t-tau1)*Ac(t-tau1) - k2*A*R + k3*AR',
    'R'  : 'qR - k2*A*R + k3*AR',
    'AR' : 'k2*A*R - k3*AR'
    }

# define parameters
params = {
    'k1' : k1_val,           
    'k2' : 0.5,
    'k3' : 0.001,
    'qR' : qR_val,
    'tau1': 40. 
    }

# initial conditions
init_cond = {
    'Aci' : 500,
    'Ai' : 0.,
    'Ri' : 0.,
    'ARi' : 0.
    }

# intialize the solver
dde = dde23(eqns=eqns, params=params)

# set the simulation parameters
# (solve from t=0 to tfinal and limit the maximum step size to 1.0) 
dde.set_sim_params(tfinal = 5000, dtmax = 1.0 )

# set the history of the proteins
histfunc = {
    'Ac' : lambda t: init_cond['Aci'],
    'A' : lambda t: init_cond['Ai'], 
    'R' : lambda t: init_cond['Ri'],
    'AR' : lambda t: init_cond['ARi']    
    }
dde.hist_from_funcs(histfunc, 500)

# run the simulator
dde.run()

#sol1 = dde.sample(0,6000,0.1)


t = dde.sol['t']
Ac = dde.sol['Ac']
A = dde.sol['A']
R = dde.sol['R']
AR = dde.sol['AR']

# print the IC's and constants
print("The parameters used were:")
print(json.dumps(params, indent = 1, sort_keys = True))

print("The initial conditions used were: ")
print(json.dumps(init_cond, indent = 1, sort_keys = True))
print 'k1*tau*e = ', params['tau1']*params['k1']*np.e


plt.figure(1)
plt.title("Reg aPKC DDE System")
pl.xlim([t[0],t[-1]])
plt.plot(t, R,'r', label = "$Reg$")
plt.plot(t,A, 'g', label = "$aPKC$")
plt.ylabel("$\mu$M")
plt.xlabel("time (s)")
pl.legend()


plt.figure(2)
plt.title("AR")
pl.xlim([t[0],t[-1]])
plt.plot(t, AR)
plt.ylabel("$\mu$M")
plt.xlabel("time (s)")
pl.legend()

plt.figure(3)
pl.xlim([t[0],t[-1]])
pl.plot(t,k1_val*R*Ac)
pl.axhline(y=qR_val, hold=None, c="red")
pl.legend()

plt.show()


