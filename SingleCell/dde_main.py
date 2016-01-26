import pdb

# biochemical signaling model 

# import pydelay and numpy and pylab
import numpy as np
from scipy.integrate import odeint
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
    dde.set_sim_params(tfinal=tf, dtmax=1)

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


def fun(y, t, Reg,tspan):
    
    # set tolerance to find indice of time to use in Reg
    tol = 25 
    
    # Find index to be used in Reg function
    index = np.where(np.logical_and(tspan >= t - tol, tspan <= t + tol))[0][0]

    Reg_cur = Reg[index]

    f = const.k_plus * Reg_cur  
    
    return f

def myosin_conc(y0,Reg,t):
    tspan = t
    return  odeint(fun,y0,t,args=(Reg,tspan,))




