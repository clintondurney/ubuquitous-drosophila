# biochemical signaling model 

from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from pydelay import dde23
import csv
import json

# the model equations 
eqns = { 
    'Ac' : '-k1*Rm*Ac',
    'Am' : 'k1*Rm(t-tau1)*Ac(t-tau1) - k2*Am*Rm + k3*AR - k2*Am*Bm + k3*AB',
    'Bc' : '-k4*Am*Bc',
    'Bm' : 'k4*Am(t-tau2)*Bc(t-tau2) - k2*Am*Bm + k3*AB',
    'Rm' : 'qR - k2*Am*Rm + k3*AR',
    'AB' : 'k2*Am*Bm - k3*AB',
    'AR' : 'k2*Am*Rm - k3*AR'
    }

# define parameters
params = {
    'k1' : 0.1,           
    'k2' : 30.75,
    'k3' : 0.0001,
    'k4' : 0.1,
    'qR' : 0.0004,
    'tau1' : 48.0,
    'tau2' : 59.3
    }

# initial conditions
init_cond = {
    'Ac' : 4.,
    'Am' : 0,
    'Bc' : 1.,
    'Bm' : 0,
    'Rm' : 0.,
    'AB' : 0,
    'AR' : 0
    }

# intialize the solver
dde = dde23(eqns=eqns, params=params)

# set the simulation parameters
# (solve from t=0 to tfinal and limit the maximum step size to 1.0) 
dde.set_sim_params(tfinal=10000, dtmax = 1.0 )


# set the history of the proteins
histfunc = {
    'Ac' : lambda t: init_cond['Ac'], 
    'Am' : lambda t: init_cond['Am'],
    'Bc' : lambda t: init_cond['Bc'],
    'Bm' : lambda t: init_cond['Bm'],
    'Rm' : lambda t: init_cond['Rm'],
    'AB' : lambda t: init_cond['AB'],
    'AR' : lambda t: init_cond['AR']
        }
dde.hist_from_funcs(histfunc, 500)

# run the simulator
dde.run()

sol1 = dde.sample(0,10000,0.1)

t = sol1['t']
Ac = sol1['Ac']
Am = sol1['Am']
Bc = sol1['Bc']
Bm = sol1['Bm']
Rm = sol1['Rm']
AB = sol1['AB']
AR = sol1['AR']

# print the IC's and constants
print("The parameters used were:")
print(json.dumps(params, indent = 1, sort_keys=True))

print("The initial conditions used were: ")
print(json.dumps(init_cond, indent = 1, sort_keys=True))

# Create plots 
plt.figure(1)
plt.title("$Reg_m$")
plt.plot(t, Rm,'r')
pl.xlim([t[0],t[-1]])
plt.ylim(0,0.1)
plt.xlabel("Time (s)")
plt.ylabel("Concentration ($\mu M$)")

plt.figure(2)
plt.title("$aPKC_m$")
plt.plot(t,Am,'g')
pl.xlim([t[0],t[-1]])
pl.ylim([0,0.1])
plt.xlabel("Time (s)")
plt.ylabel("Concentration ($\mu M$)")

plt.figure(3)
plt.title("$Baz_m$")
plt.plot(t,Bm,'b')
pl.xlim([t[0],t[-1]])
pl.ylim([0,0.1])
plt.xlabel("Time (s)")
plt.ylabel("Concentration ($\mu M$)")

plt.figure(4)
plt.title("Medial Reg, aPKC, Baz")
pl.xlim([t[0],t[-1]])
pl.ylim([0,1])
plt.plot(t, Rm/np.amax(Rm),'r', label = "$Reg_m$")
plt.plot(t,Am/np.amax(Am), 'g', label = "$aPKC_m$")
plt.plot(t,Bm/np.amax(Bm), 'b', label = "$Baz_m$")
plt.ylabel("Normalized Concentration")
plt.xlabel("time (s)")
pl.legend()

plt.show()

out = np.column_stack((t,Rm,Am,Bm))
np.savetxt("output.csv",out,delimiter=",",header='t,Rm,Am,Bm',comments='')


