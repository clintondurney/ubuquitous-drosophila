import pdb
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import globals as const
from dde_main import *
from scipy.spatial import distance
import csv

###########
#
# main4.py
#
#
# Author: Clinton H. Durney
# Email: cdurney@math.ubc.ca
#
# Last Edit: 5/14/16
#
# To Do:    
#
#
###########

# Set-up output file
BioParamsFile = open('BioParams.csv','w')
BioParamsWriter = csv.writer(BioParamsFile,delimiter='\t')
BioParamsWriter.writerow(["time", "Reg","myo0_1","myo1_17","myo1_7","myo17_13","myo1_10","myo10_17"])

G = nx.Graph()
H = nx.Graph()

r = const.l0

# Left most cell set up
origin = (0.0,0.0)
p0 = (r,0.0)
p1 = (r*np.cos(np.pi/3),r*np.sin(np.pi/3))
p2 = (r*np.cos(2*np.pi/3),r*np.sin(2*np.pi/3))
p3 = (r*np.cos(3*np.pi/3),r*np.sin(3*np.pi/3))
p4 = (r*np.cos(4*np.pi/3),r*np.sin(4*np.pi/3))
p5 = (r*np.cos(5*np.pi/3),r*np.sin(5*np.pi/3))
nodes = [p0,p1,p2,p3,p4,p5]
   
nodes = [p0,p1,p2,p3,p4,p5]
i = 0
G.add_node(i,pos=origin, center=True, phase_angle=0, boundary=[1,2,3,4,5,6])
i += 1 
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# Middle Top Cell with Cap
origin = 1.5*r, r*np.cos(np.pi/6)
p1 = (origin[0] + r*np.cos(np.pi/3), origin[1] +r* np.sin(np.pi/3))
p2 = (origin[0] + r*np.cos(2*np.pi/3), origin[1] + r*np.sin(2*np.pi/3))
    
nodes = [p1,p2]
    
G.add_node(i,pos=origin,center=True, phase_angle = 10000, boundary=[16,8,9,2,1,17])
i +=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# Middle Bottom Cell with Cap
origin = 1.5*r, -r*np.cos(np.pi/6)
p4 = (origin[0] + r*np.cos(4*np.pi/3), origin[1] + r*np.sin(4*np.pi/3))
p5 = (origin[0] + r*np.cos(5*np.pi/3), origin[1] + r*np.sin(5*np.pi/3))
nodes = [p4,p5]
    
G.add_node(i,pos=origin,center=True,phase_angle = 10000, boundary=[18,17,1,6,11,12])
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
    
G.add_node(i,pos=origin,center=True,phase_angle = 0, boundary=[14,15,16,17,18,19])
i+=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# edges that are only passively elastic
G.add_path([1,2,3,4,5,6,1],beta=0,myosin=0, color = 'r')
G.add_path([1,17,16,8,9,2],beta=0,myosin=0, color = 'r')
G.add_path([6,11,12,18,17,1],beta=0,myosin=0,color = 'r')
G.add_path([18,19,14,15,16],beta=0,myosin=0, color = 'r')
    
# edges that have active force component
G.add_edges_from([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(7,16),(7,8),(7,9),(7,2),(7,1),(7,17)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(10,18),(10,17),(10,1),(10,6),(10,11),(10,12)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(13,14),(13,15),(13,16),(13,17),(13,18),(13,19)],beta=const.beta,myosin=const.myo0)

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

for index in range(160000,len(dt)):
    if t[index] >= 0:
        H = G.copy() 
        ## Update myosin concentration on each spoke ##
        for center in G.nodes_iter(data=True):
            if center[1]['center']==True:
                # Calculate cell area
                outer = [G.node[element]['pos'] for element in center[1]['boundary']]
                cell_area = CellArea(outer)
                
                for neighbor in G.neighbors(center[0]):
                    # Calculate area of adjacent triangles of the spoke    
                    inner = [G.node[neighbor]['pos']]
                    temp = list(set(G.neighbors(center[0])) & set(G.neighbors(neighbor)))
                    inner.append(G.node[temp[0]]['pos'])
                    inner.append(center[1]['pos'])
                    inner.append(G.node[temp[1]]['pos'])
                    spoke_area = CellArea(inner)
                    geo_frac = spoke_area/cell_area
                    
                    # Calculate necessary parameters for dm/dt
                    length = distance.euclidean(G.node[center[0]]['pos'],G.node[neighbor]['pos'])
                    myosin_current = G[center[0]][neighbor]['myosin']
                    
                    #update myosin on this edge
                    G[center[0]][neighbor]['myosin'] = dmyosin(myosin_current, geo_frac*Rm[index-center[1]['phase_angle']], length, dt[index])


        ## Update force ##
        # iterate over all nodes in graph
        for point in G.nodes_iter():
    	    # iterate over all neighbors of node
    	    total_force = [0,0]
    	    for neighbor in G.neighbors(point):
                # calculate the unit vector from node to neighbor
                dir_vector = unit_vector(H.node[point]['pos'],H.node[neighbor]['pos'])
                # calculate magnitude of force
                length = distance.euclidean(H.node[point]['pos'],H.node[neighbor]['pos'])
                mag_force = calc_force(length,G[point][neighbor]['myosin'])
                total_force = np.sum([total_force,mag_force*np.array(dir_vector)],axis=0)
            
            if point not in [2,3,4,5,6,11,12,18,19,14,15,16,8,9]:
                G.node[point]['pos'] = d_pos(H.node[point]['pos'],total_force, dt[index])

        if index % 100 == 0:
            print t[index]
            history.append(G.copy())
            
            BioParamsWriter.writerow([t[index], Rm[index], G[0][1]['myosin'], G[1][17]['myosin'], G[1][7]['myosin'], G[17][13]['myosin'],
                G[1][10]['myosin'],G[10][17]['myosin']])


        if index % 3000  == 0:
            for i in range(0,len(history)):
                plt.clf()
                pos = nx.get_node_attributes(history[i],'pos')

                edges,colors = zip(*nx.get_edge_attributes(history[i],'color').items())
                nx.draw(history[i],pos,edgelist=edges,edge_color=colors,width=2)

                plt.xlim(-10,25)
                plt.ylim(-12,12)
                plt.axis("on")
                plt.grid("on")
                
                # plt.show()
                plt.savefig('tmp%03d.png'%i)





