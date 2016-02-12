import pdb

# biochemical signaling model 

# import pydelay and numpy and pylab
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import pylab as pl
from pydelay import dde23
import globals as const
import csv
import json

def dde_initializer(Ac_i,Am_i,Bc_i,Bm_i,Rm_i,AB_i,AR_i,tf):
    # the model equations 
    eqns = { 
        'Ac' : '-k1*Rm*Ac',
        'Am' : 'k1*Rm(t-tau1)*Ac(t-tau1)*Heavi(t-tau1) - k2*Am*Rm + k2*AR - k2*Am*Bm + k3*AB',
        'Bc' : '-k4*Am*Bc',
        'Bm' : 'k4*Am(t-tau2)*Bc(t-tau2)*Heavi(t-tau2) - k2*Am*Bm + k3*AB',
        'Rm' : 'qR - k2*Am*Rm + k2*AR',
        'AB' : 'k2*Am*Bm - k3*AB',
        'AR' : 'k2*Am*Rm - k2*AR',
        }

    # define parameters
    params = {
        'k1' : const.k1,           
        'k2' : const.k2,
        'k3' : const.k3,
        'k4' : const.k4,
        'qR' : const.qR,
        'tau1' : const.tau1,
        'tau2' : const.tau2,
        }

    # initial conditions
    init_cond = {
        'Ac' : Ac_i,
        'Am' : Am_i,
        'Bc' : Bc_i,
        'Bm' : Bm_i,
        'Rm' : Rm_i,
        'AB' : AB_i,
        'AR' : AR_i,
        }
    
    # intialize the solver
    dde = dde23(eqns=eqns, params=params)

    # set the simulation parameters
    # (solve from t=0 to t=tf and limit the maximum step size to 1.0) 
    dde.set_sim_params(tfinal=tf, dtmax=.01)

    # set the history of the proteins
    histfunc = {
        'Ac' : lambda t: init_cond['Ac'], 
        'Am' : lambda t: init_cond['Am'],
        'Bc' : lambda t: init_cond['Bc'],
        'Bm' : lambda t: init_cond['Bm'],
        'Rm' : lambda t: init_cond['Rm'],
        'AB' : lambda t: init_cond['AB'],
        'AR' : lambda t: init_cond['AR'],
        }
    
    dde.hist_from_funcs(histfunc,500)

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

    return(t,Ac,Am,Bc,Bm,Rm,AB,AR) 


def dde_solver(t_i,Ac_i,Am_i,Bc_i,Bm_i,Rm_i,AB_i,AR_i,tf):
    # the model equations 
    eqns = { 
        'Ac' : '-k1*Rm*Ac',
        'Am' : 'k1*Rm(t-tau1)*Ac(t-tau1)*Heavi(t-tau1) - k2*Am*Rm + k2*AR - k2*Am*Bm + k3*AB',
        'Bc' : '-k4*Am*Bc',
        'Bm' : 'k4*Am(t-tau2)*Bc(t-tau2)*Heavi(t-tau2) - k2*Am*Bm + k3*AB',
        'Rm' : 'qR - k2*Am*Rm + k2*AR',
        'AB' : 'k2*Am*Bm - k3*AB',
        'AR' : 'k2*Am*Rm - k2*AR',
        }

    # define parameters
    params = {
        'k1' : const.k1,           
        'k2' : const.k2,
        'k3' : const.k3,
        'k4' : const.k4,
        'qR' : const.qR,
        'tau1' : const.tau1,
        'tau2' : const.tau2,
        }

    # initial conditions
    histdict = {
        't'  : t_i,
        'Ac' : Ac_i,
        'Am' : Am_i,
        'Bc' : Bc_i,
        'Bm' : Bm_i,
        'Rm' : Rm_i,
        'AB' : AB_i,
        'AR' : AR_i,
        }
    
    # intialize the solver
    dde = dde23(eqns=eqns, params=params)

    dde.hist_from_arrays(histdict,useend=False)
    # set the simulation parameters
    
    # (solve from t=0 to t=tf and limit the maximum step size to 1.0) 
    dde.set_sim_params(tfinal=tf, dtmax=1)

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

    return(t,Ac,Am,Bc,Bm,Rm,AB,AR)

###############
def normed_direction(A,B):
    # Calculate the normed vector BA

    dist = distance.euclidean(A,B)

    return ((A[0]-B[0])/dist,(A[1]-B[1])/dist)

###############
def d_pos(x_old,y_old,f,direction,dt):
    
    x_new = x_old + (dt/const.eta)*f*direction[0]

    y_new = y_old + (dt/const.eta)*f*direction[1]

    return [x_new,y_new]


###############
def calc_force(l, myosin):
    
    return const.mu*(l-const.l0) + const.beta*myosin

###############
def kminus(myo,length):
    
    return(const.k1*np.exp(-const.k2*(calc_force(length,myo))))

################
def fun(y,signal,length):
    
    geo_frac = 1                # evenly distributed to start (single cell)
    
    return(const.k_plus*signal*geo_frac - kminus(y,length)*y)

################
def dmyosin(y,Reg,length, dt):

    return(y+dt*fun(y,Reg,length))




