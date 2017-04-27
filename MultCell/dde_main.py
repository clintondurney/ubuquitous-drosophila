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
    def gen_nodes(ori):
        nodes = []
        for n in range(0,6):
            nodes.append((ori[0] + r*np.cos(n*np.pi/3), ori[1] + r*np.sin(n*np.pi/3)))
        return nodes

    def add_nodes(nodes, i):
        pos = nx.get_node_attributes(G,'pos')
        cen_index = i-1
        centers.append(cen_index)
        AS_boundary = []
        spokes = []
        for node in nodes:
            add_node = True
            for existing_node in range(0,len(pos)):
                if distance.euclidean(pos[existing_node],node) < 10**(-7):
                    add_node = False
                    AS_boundary.append(existing_node)
                    spokes.append((cen_index,existing_node))
                    break

            if add_node == True:
                G.add_node(i,pos=node)
                i += 1
                AS_boundary.append(i-1)
                spokes.append((cen_index,i-1))

        return AS_boundary, spokes, i

    def add_spokes_edges(spokes, boundary):
        boundary.append(boundary[0])
        G.add_edges_from(spokes,beta=10,myosin=1000)    
        G.add_path(boundary,beta=0,myosin=0,color='r')

        return

    G = nx.Graph()

    r = const.l_initial        # initial spoke length
    num_cells = const.num_center_row  # number of cells in center row

    centers = []
    i = 0
    # Center cell set up
    origin = (0.0,0.0)
    G.add_node(i,pos=origin)
    i += 1

    nodes = gen_nodes(origin)
    AS_boundary, spokes, i = add_nodes(nodes,i)
    add_spokes_edges(spokes, AS_boundary)

    for index in range(1,int((num_cells - 1)/2.)+1):
        # # Step Up
        origin = (0, np.sqrt(3)*r*index)
        G.add_node(i,pos=origin)
        i += 1

        nodes = gen_nodes(origin)
        AS_boundary, spokes, i = add_nodes(nodes,i)
        add_spokes_edges(spokes, AS_boundary)

        # # # Step down
        origin = (0, -np.sqrt(3)*r*index)
        G.add_node(i,pos=origin)
        i += 1

        nodes = gen_nodes(origin)
        AS_boundary, spokes, i = add_nodes(nodes,i)
        add_spokes_edges(spokes, AS_boundary)

    for index in range(1,num_cells):  
        if (num_cells - index) % 2 == 0:
            for j in range(1,(num_cells-index),2):
                origin = ((3/2.)*r*index,(np.sqrt(3)/2.)*r*j)
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)

                origin = ((3/2.)*r*index,(-np.sqrt(3)/2.)*r*j)
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)

            # Step Left

                origin = (-(3/2.)*r*index,(np.sqrt(3)/2.)*r*j)
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)

                origin = (-(3/2.)*r*index,(-np.sqrt(3)/2.)*r*j)
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)

        else:
            for j in range(0,(num_cells-index),2):
                origin = (3*(1/2.)*r*index, (np.sqrt(3)/2.)*r*j)
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)
                
                if j != 0:
                    origin = (3*(1/2.)*r*index, -(np.sqrt(3)/2.)*r*j)
                    G.add_node(i,pos=origin)
                    i += 1

                    nodes = gen_nodes(origin)
                    AS_boundary, spokes, i = add_nodes(nodes,i)
                    add_spokes_edges(spokes, AS_boundary)

                # Step Left
                origin = (-3*(1/2.)*r*index, (np.sqrt(3)/2.)*r*j)
                G.add_node(i,pos=origin)
                i += 1

                nodes = gen_nodes(origin)
                AS_boundary, spokes, i = add_nodes(nodes,i)
                add_spokes_edges(spokes, AS_boundary)
                
                if j != 0:
                    origin = (-3*(1/2.)*r*index, -(np.sqrt(3)/2.)*r*j)
                    G.add_node(i,pos=origin)
                    i += 1

                    nodes = gen_nodes(origin)
                    AS_boundary, spokes, i = add_nodes(nodes,i)
                    add_spokes_edges(spokes, AS_boundary)

    nx.set_node_attributes(G, 'time_lag', 0)
    
    for j in centers:
            G.node[j]['time_lag'] = np.random.randint(0,2000)
    
    AS_boundary = []
    for j in G.nodes_iter():
        if G.degree(j) == 3 or G.degree(j) == 5:
            AS_boundary.append(j)
    
    temp_sort = [(AS_boundary[0],0)]
    pos = nx.get_node_attributes(G,'pos')
    u = unit_vector(pos[0],pos[AS_boundary[0]])
    for index in range(1,len(AS_boundary)):
            v = unit_vector(pos[0],pos[AS_boundary[index]])
            dot = np.dot(u,v)
            det = np.linalg.det([u,v])
            angle = np.arctan2(det,dot)     
            temp_sort.append((AS_boundary[index],angle))
    temp_sort = sorted(temp_sort, key=lambda tup: tup[1])
    AS_boundary = [temp_sort[j][0] for j in range(0,len(temp_sort))]
    
    epidermis = []

    for index in range(0,len(AS_boundary)):
            temp = list(set(centers).intersection(G.neighbors(AS_boundary[index]))) 
            if len(temp) == 1:
                dirn = unit_vector(pos[temp[0]],pos[AS_boundary[index]])
            else:
                v_1 = unit_vector(pos[AS_boundary[index]],pos[temp[0]])
                v_2 = unit_vector(pos[AS_boundary[index]],pos[temp[1]])
                dirn = -(v_1[0]+v_2[0]),-(v_1[1]+v_2[1])
            x = pos[AS_boundary[index]][0] + 10*dirn[0]
            y = pos[AS_boundary[index]][1] + 10*dirn[1]
            G.add_node(i,pos=(x,y))
            G.add_edges_from([(AS_boundary[index],i)],beta=0,myosin=0, color='b')
            epidermis.append(i)
            i += 1
    epidermis.append(epidermis[0]) 
    G.add_path(epidermis,beta=0,myosin=0,color='b')
    
    return G, centers, epidermis
