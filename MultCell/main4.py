import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import globals as const
from dde_main import *
from scipy.spatial import distance


###########
#
# main4.py
#
#
# Author: Clinton H. Durney
# Email: cdurney@math.ubc.ca
#
# Last Edit: 04/28/16
#
# To Do:    
#
#
###########

G = nx.Graph()

r = const.initialspokelength

# Left most cell set up
origin = (0,0)
p0 = (r,0)
p1 = (r*np.cos(np.pi/3),r*np.sin(np.pi/3))
p2 = (r*np.cos(2*np.pi/3),r*np.sin(2*np.pi/3))
p3 = (r*np.cos(3*np.pi/3),r*np.sin(3*np.pi/3))
p4 = (r*np.cos(4*np.pi/3),r*np.sin(4*np.pi/3))
p5 = (r*np.cos(5*np.pi/3),r*np.sin(5*np.pi/3))
nodes = [p0,p1,p2,p3,p4,p5]
   
nodes = [p0,p1,p2,p3,p4,p5]
i = 0
G.add_node(i,pos=origin, center=True)
i += 1 
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# Middle Top Cell with Cap
origin = 1.5*r, r*np.cos(np.pi/6)
p1 = (origin[0] + r*np.cos(np.pi/3), origin[1] +r* np.sin(np.pi/3))
p2 = (origin[0] + r*np.cos(2*np.pi/3), origin[1] + r*np.sin(2*np.pi/3))
    
nodes = [p1,p2]
    
G.add_node(i,pos=origin,center=True)
i +=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# Middle Bottom Cell with Cap
origin = 1.5*r, -r*np.cos(np.pi/6)
p4 = (origin[0] + r*np.cos(4*np.pi/3), origin[1] + r*np.sin(4*np.pi/3))
p5 = (origin[0] + r*np.cos(5*np.pi/3), origin[1] + r*np.sin(5*np.pi/3))
nodes = [p4,p5]
    
G.add_node(i,pos=origin,center=True)
i+=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# Right most cell
origin = 1.5*2*r, 0
p0 = (origin[0] + r, origin[1] + 0)
p1 = (origin[0] + r*np.cos(np.pi/3), origin[1] + r*np.sin(np.pi/3))
p2 = (origin[0] + r*np.cos(2*np.pi/3), origin[1] + r*np.sin(2*np.pi/3))
p3 = (origin[0] + r*np.cos(3*np.pi/3), origin[1] + r*np.sin(3*np.pi/3))
p4 = (origin[0] + r*np.cos(4*np.pi/3), origin[1] + r*np.sin(4*np.pi/3))
p5 = (origin[0] + r*np.cos(5*np.pi/3), origin[1] + r*np.sin(5*np.pi/3))
nodes = [p0,p1,p2,p3,p4,p5]
    
G.add_node(i,pos=origin,center=True)
i+=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# edges that are only passively elastic
G.add_path([1,2,3,4,5,6,1],beta=0)
G.add_path([1,17,16,8,9,2],beta=0)
G.add_path([6,11,12,18,17,1],beta=0)
G.add_path([18,19,14,15,16],beta=0)
    
# edges that have active force component
G.add_edges_from([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6)],beta=10,myosin=const.myo0)
G.add_edges_from([(7,16),(7,8),(7,9),(7,2),(7,1),(7,17)],beta=10,myosin=const.myo0)
G.add_edges_from([(10,18),(10,17),(10,1),(10,6),(10,11),(10,12)],beta=10,myosin=const.myo0)
G.add_edges_from([(13,14),(13,15),(13,16),(13,17),(13,18),(13,19)],beta=10,myosin=const.myo0)
    

history = [G.copy()]

# Initial conditions of Biochemical parameters
Ac = const.Ac0
Am = const.Am0
Bc = const.Bc0
Bm = const.Bm0
Rm = const.Rm0
AB = const.AB0
AR = const.AR0

# Run dde solver for the Biochemical concentrations
tf = 6000
(t,Ac,Am,Bc,Bm,Rm,AB,AR) = dde_initializer(Ac,Am,Bc,Bm,Rm,AB,AR,tf)

dt = np.diff(t)

for index in range(1,len(dt)):
    if t[index] >= 0:
        
        # Update myosin concentration on each spoke
        for center in G.nodes_iter(data=True):
            if center[1]['center']==True:
                for neighbor in G.neighbors(center[0]):
                    length = distance.euclidean(G.node[center[0]]['pos'],G.node[neighbor]['pos'])
                    myosin_current = G[center[0]][neighbor]['myosin']
                    G[center[0]][neighbor]['myosin'] = dmyosin(myosin_current, Rm[index], length, dt[index])


        # Update force
        for point in G.nodes_iter(data=True):

    
        # Update location

        




