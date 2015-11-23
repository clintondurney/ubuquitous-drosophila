# biochemical signaling model 

# import pydelay and numpy and pylab
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from pydelay import dde23
import csv
import json

# the model equations 
eqns = { 
    'Ac' : '-k1*Rm*Ac',
    'Am' : 'k1*Rm(t-tau1)*Ac(t-tau1)*Heavi(t-tau1) - k2*Am*Rm + k2*AR - k2*Am*Bm + k3*AB',
    'Bc' : '-k4*Am*Bc',
    'Bm' : 'k4*Am(t-tau2)*Bc(t-tau2)*Heavi(t-tau2) - k2*Am*Bm + k3*AB',
    'Rm' : 'qR - k2*Am*Rm + k2*AR',
    'AB' : 'k2*Am*Bm - k3*AB',
    'AR' : 'k2*Am*Rm - k2*AR'
    }

# define parameters
params = {
    'k1' : 0.000003,           
    'k2' : 0.15,
    'k3' : 0.009,
    'k4' : 0.000003,
    'qR' : 8,
    'tau1' : 75,
    'tau2' : 80
    }

# initial conditions
init_cond = {
    'Ac' : 60000.,
    'Am' : 0,
    'Bc' : 4000.,
    'Bm' : 0,
    'Rm' : 500.,
    'AB' : 0,
    'AR' : 0
    }

# intialize the solver
dde = dde23(eqns=eqns, params=params)

# set the simulation parameters
# (solve from t=0 to tfinal and limit the maximum step size to 1.0) 
dde.set_sim_params(tfinal=6000, dtmax=0.1)

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

# get the solutions from the history dict
t = dde.sol['t']
Ac= dde.sol['Ac']
Am= dde.sol['Am']
Bc = dde.sol['Bc']
Bm = dde.sol['Bm']
Rm = dde.sol['Rm']
AB = dde.sol['AB']
AR = dde.sol['AR']

# print the IC's and constants
print "The parameters used were: "
print json.dumps(params, indent = 1)

print "The initial conditions used were: "
print json.dumps(init_cond, indent = 1)

# Create plots 
plt.figure(1)
plt.title("$Reg_m$")
plt.plot(t, Rm,'r')
pl.xlim([t[0],t[-1]])
plt.xlabel("Time (s)")
plt.ylabel("$Reg_m$")

plt.figure(2)
plt.title("$aPKC_m$")
plt.plot(t,Am,'g')
pl.xlim([t[0],t[-1]])
plt.xlabel("Time (s)")
plt.ylabel("$aPKC_m$")

plt.figure(3)
plt.title("$Baz_m$")
plt.plot(t,Bm,'b')
pl.xlim([t[0],t[-1]])
plt.xlabel("Time (s)")
plt.ylabel("$Baz_m$")

plt.figure(4)
plt.title("Medial Reg, aPKC and Baz vs. Time")
pl.xlim([t[0],t[-1]])
plt.plot(t, Rm/np.amax(Rm),'r', label = "$Reg_m$")
plt.plot(t,Am/np.amax(Am), 'g', label = "$aPKC_m$")
plt.plot(t,Bm/np.amax(Bm), 'b', label = "$Baz_m$")
pl.legend()

plt.show()




