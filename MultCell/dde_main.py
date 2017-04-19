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
import networkx as nx
import itertools



def dde_initializer(Ac_i,Am_i,Bc_i,Bm_i,Rm_i,AB_i,AR_i,tf,dt):
    # the model equations 
    eqns = { 
        'Ac' : '-k1*Rm*Ac',
        'Am' : 'k1*Rm(t-tau1)*Ac(t-tau1)*Heavi(t-tau1) - k2*Am*Rm + k3*AR - k2*Am*Bm + k3*AB',
        'Bc' : '-k4*Am*Bc',
        'Bm' : 'k4*Am(t-tau2)*Bc(t-tau2)*Heavi(t-tau2) - k2*Am*Bm + k3*AB',
        'Rm' : 'qR - k2*Am*Rm + k3*AR',
        'AB' : 'k2*Am*Bm - k3*AB',
        'AR' : 'k2*Am*Rm - k3*AR',
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
    dde.set_sim_params(tfinal=tf, dtmax=1.0)

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

    # Get solution at every t=0.1
    sol1 = dde.sample(0,tf,0.1)

    # get the solutions from the history dict
#    t = dde.sol['t']
#    Ac= dde.sol['Ac']
#    Am= dde.sol['Am']
#    Bc = dde.sol['Bc']
#    Bm = dde.sol['Bm']
#    Rm = dde.sol['Rm']
#    AB = dde.sol['AB']
#    AR = dde.sol['AR']

    t = sol1['t']
    Ac= sol1['Ac']
    Am= sol1['Am']
    Bc = sol1['Bc']
    Bm = sol1['Bm']
    Rm = sol1['Rm']
    AB = sol1['AB']
    AR = sol1['AR']

    return(t,Ac,Am,Bc,Bm,Rm,AB,AR) 


def dde_solver(t_i,Ac_i,Am_i,Bc_i,Bm_i,Rm_i,AB_i,AR_i,tf):
    # the model equations 
    eqns = { 
        'Ac' : '-k1*Rm*Ac',
        'Am' : 'k1*Rm(t-tau1)*Ac(t-tau1)*Heavi(t-tau1) - k2*Am*Rm + k3*AR - k2*Am*Bm + k3*AB',
        'Bc' : '-k4*Am*Bc',
        'Bm' : 'k4*Am(t-tau2)*Bc(t-tau2)*Heavi(t-tau2) - k2*Am*Bm + k3*AB',
        'Rm' : 'qR - k2*Am*Rm + k3*AR',
        'AB' : 'k2*Am*Bm - k3*AB',
        'AR' : 'k2*Am*Rm - k3*AR',
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
def unit_vector(A,B):
    # Calculate the unit vector from A to B 

    dist = distance.euclidean(A,B)

    return ((B[0]-A[0])/dist,(B[1]-A[1])/dist)

###############
def d_pos(position,force,dt):
    
    x_new = position[0] + (dt/const.eta)*force[0]

    y_new = position[1] + (dt/const.eta)*force[1]

    return (x_new,y_new)


###############
def calc_force(l, myosin):
    
    return const.mu*(l-const.l0) + const.beta*myosin

###############
def kminus(myo,length):
    
    if (calc_force(length,myo)) > 0:
        return(const.c_1*np.exp(-const.c_2*(calc_force(length,myo))))
    else:
        return const.c_1

################
def fun(y,signal,length):
    # needs updated for geo distribution
    
    return(const.k_plus*signal - kminus(y,length)*y)

################
def dmyosin(y,Reg,length, dt):

    return(y+dt*fun(y,Reg,length))

################
def CellArea(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

################
def conc2mol(IC, cell_vol):
    Mol = IC*10**(-3)*10**(-12)*6.022*10**(23)*cell_vol
    return Mol

################
def tissue():
    r = const.l_initial
    xx = np.arange(0,15*r,r)
    yy = np.arange(0,11*r*np.sqrt(3)*0.5,r*np.sqrt(3)*0.5)
    points = list(itertools.product(xx,yy))

    # shear the points
    m = np.sqrt(3)/3

    new_points = []
    for i in range(0,len(points)):
        x = points[i][0] + m*points[i][1]
        y = points[i][1]
        new_points.append((x,y))  

    points = new_points

    G = nx.Graph()
    for i in range(0,len(points)):
        G.add_node(i,pos=points[i])

    section_1 = [0,1,2,3,4,11,12,13,14,22,23,24,25,33,34,35,44,45,46,55,56,66,67,77,88]
    section_2 = [121,132,133,143,144,145,154,155,156,157]
    section_3 = [7,8,9,10,19,20,21,31,32,43]
    section_4 = [76,87,98,109,120,131,142,153,164,97,108,119,130,141,152,163,118,129,140,151,162,139,150,161,160]

    G.remove_nodes_from(section_1)
    G.remove_nodes_from(section_2)
    G.remove_nodes_from(section_3)
    G.remove_nodes_from(section_4)

    G = nx.convert_node_labels_to_integers(G,first_label=0)

    G.add_path([0,1,4,7,6,2,0],beta=0,myosin=0,color='r')
    G.add_path([6,11,18,19,13,7],beta=0,myosin=0,color='r')
    G.add_path([18,26,35,36,28,19],beta=0,myosin=0,color='r')
    G.add_path([4,5,9,14,13],beta=0,myosin=0,color='r')
    G.add_path([14,21,29,28],beta=0,myosin=0,color='r')
    G.add_path([9,10,16,22,21],beta=0,myosin=0,color='r')
    G.add_path([22,31,39,38,29],beta=0,myosin=0,color='r')
    G.add_path([38,46,45,36],beta=0,myosin=0,color='r')
    G.add_path([35,43,52,53,45],beta=0,myosin=0,color='r')
    G.add_path([31,32,24,17,16],beta=0,myosin=0,color='r')
    G.add_path([32,41,42,34,25,24],beta=0,myosin=0,color='r')
    G.add_path([39,48,49,41],beta=0,myosin=0,color='r')
    G.add_path([46,55,56,48],beta=0,myosin=0,color='r')
    G.add_path([53,62,63,55],beta=0,myosin=0,color='r')
    G.add_path([52,60,69,70,62],beta=0,myosin=0,color='r')
    G.add_path([70,77,78,72,63],beta=0,myosin=0,color='r')
    G.add_path([72,73,65,56],beta=0,myosin=0,color='r')
    G.add_path([65,66,58,49],beta=0,myosin=0,color='r')
    G.add_path([58,59,51,42],beta=0,myosin=0,color='r')
    G.add_path([66,75,76,68,59],beta=0,myosin=0,color='r')
    G.add_path([73,80,81,75],beta=0,myosin=0,color='r')
    G.add_path([78,84,85,80],beta=0,myosin=0,color='r')
    G.add_path([81,87,88,83,76],beta=0,myosin=0,color='r')
    G.add_path([85,89,90,87],beta=0,myosin=0,color='r')
    G.add_path([90,93,94,92,88],beta=0,myosin=0,color='r')


    G.add_edges_from([(3,0),(3,1),(3,4),(3,7),(3,6),(3,2)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(8,14),(8,9),(8,5),(8,4),(8,7),(8,13)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(12,19),(12,13),(12,7),(12,6),(12,11),(12,18)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(15,22),(15,16),(15,10),(15,9),(15,14),(15,21)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(20,29),(20,21),(20,14),(20,13),(20,19),(20,28)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(27,36),(27,28),(27,19),(27,18),(27,26),(27,35)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(23,32),(23,24),(23,17),(23,16),(23,22),(23,31)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(30,39),(30,31),(30,22),(30,21),(30,29),(30,38)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(37,46),(37,38),(37,29),(37,28),(37,36),(37,45)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(44,53),(44,45),(44,36),(44,35),(44,43),(44,52)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(33,42),(33,34),(33,25),(33,24),(33,32),(33,41)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(40,49),(40,41),(40,32),(40,31),(40,39),(40,48)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(47,56),(47,48),(47,39),(47,38),(47,46),(47,55)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(54,63),(54,55),(54,46),(54,45),(54,53),(54,62)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(61,70),(61,62),(61,53),(61,52),(61,60),(61,69)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(50,59),(50,51),(50,42),(50,41),(50,49),(50,58)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(57,66),(57,58),(57,49),(57,48),(57,56),(57,65)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(64,73),(64,65),(64,56),(64,55),(64,63),(64,72)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(71,78),(71,72),(71,63),(71,62),(71,70),(71,77)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(67,76),(67,68),(67,59),(67,58),(67,66),(67,75)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(74,81),(74,75),(74,66),(74,65),(74,73),(74,80)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(79,85),(79,80),(79,73),(79,72),(79,78),(79,84)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(82,88),(82,83),(82,76),(82,75),(82,81),(82,87)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(86,90),(86,87),(86,81),(86,80),(86,85),(86,89)],beta=const.beta,myosin=const.myo0)
    G.add_edges_from([(91,94),(91,92),(91,88),(91,87),(91,90),(91,93)],beta=const.beta,myosin=const.myo0)

    nx.set_node_attributes(G, 'time_lag', 0)
        
    centers = [3,8,12,15,20,27,23,30,37,44,33,40,47,54,61,50,57,64,71,67,74,79,82,86,91]
                
    for j in centers:
        G.node[j]['time_lag'] = np.random.randint(0,2000)

    return G



