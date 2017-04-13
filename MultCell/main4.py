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
# Last Edit: 4/06/17
#
# To Do:    
#
#
###########

# Set-up output file
BioParamsFile = open('BioParams.csv','w')
BioParamsWriter = csv.writer(BioParamsFile,delimiter='\t')
# BioParamsWriter.writerow(["time", "Reg","myo0_1","myo1_17","myo1_7","myo17_13","myo1_10","myo10_17"])
BioParamsWriter.writerow(["time","cell_0","cell_7","cell_10","cell_13","cell_20","cell_24", "cell_28", "cell_32", "cell_35"])

G = nx.Graph()
H = nx.Graph()

r = const.l_initial

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
G.add_node(i,pos=origin, center=True, phase_angle=np.random.randint(2000), boundary=[1,2,3,4,5,6])
i += 1 
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# Middle Top Cell with Cap
origin = 1.5*r, r*np.cos(np.pi/6)
p1 = (origin[0] + r*np.cos(np.pi/3), origin[1] +r* np.sin(np.pi/3))
p2 = (origin[0] + r*np.cos(2*np.pi/3), origin[1] + r*np.sin(2*np.pi/3))
    
nodes = [p1,p2]
    
G.add_node(i,pos=origin,center=True, phase_angle = np.random.randint(2000), boundary=[16,8,9,2,1,17])
i +=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
# Middle Bottom Cell with Cap
origin = 1.5*r, -r*np.cos(np.pi/6)
p4 = (origin[0] + r*np.cos(4*np.pi/3), origin[1] + r*np.sin(4*np.pi/3))
p5 = (origin[0] + r*np.cos(5*np.pi/3), origin[1] + r*np.sin(5*np.pi/3))
nodes = [p4,p5]
    
G.add_node(i,pos=origin,center=True,phase_angle = np.random.randint(2000), boundary=[18,17,1,6,11,12])
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
    
G.add_node(i,pos=origin,center=True,phase_angle = np.random.randint(2000), boundary=[14,15,16,17,18,19])
i+=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1

###########################
origin = 1.5*2*r, 2*r*np.cos(np.pi/6)
p0 = (origin[0] + r, origin[1] + 0)
p1 = (origin[0] + r*np.cos(np.pi/3), origin[1] + r*np.sin(np.pi/3))
p2 = (origin[0] + r*np.cos(2*np.pi/3), origin[1] + r*np.sin(2*np.pi/3))
nodes = [p0,p1,p2]

G.add_node(i,pos=origin,center=True,phase_angle = np.random.randint(2000), boundary=[21,22,23,8,16,15])
i+=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1

origin = 1.5*2*r, -2*r*np.cos(np.pi/6)
p0 = (origin[0] + r, origin[1] + 0)
p4 = (origin[0] + r*np.cos(4*np.pi/3), origin[1] + r*np.sin(4*np.pi/3))
p5 = (origin[0] + r*np.cos(5*np.pi/3), origin[1] + r*np.sin(5*np.pi/3))
nodes = [p0,p4,p5]

G.add_node(i,pos=origin,center=True,phase_angle = np.random.randint(2000), boundary=[25,19,18,12,26,27])
i+=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1

origin = 3*(1.5*r), r*np.cos(np.pi/6)
p0 = (origin[0] + r, origin[1] + 0.0)
p1 = (origin[0] + r*np.cos(np.pi/3), origin[1] + r*np.sin(np.pi/3))
p5 = (origin[0] + r*np.cos(5*np.pi/3), origin[1] + r*np.sin(5*np.pi/3))
nodes = [p0,p1,p5]
    
G.add_node(i,pos=origin,center=True, phase_angle = np.random.randint(2000), boundary=[29,30,21,15,14,31])
i +=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1

origin = 3*(1.5*r), -r*np.cos(np.pi/6)
p0 = (origin[0] + r, origin[1] + 0)
p5 = (origin[0] + r*np.cos(5*np.pi/3), origin[1] + r*np.sin(5*np.pi/3))
nodes = [p0,p5]

G.add_node(i,pos=origin,center=True, phase_angle = np.random.randint(2000), boundary=[33,31,14,19,25,34])
i +=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1
    
origin = 4*(1.5*r), 0
p0 = (origin[0] + r, origin[1] + 0)
p1 = (origin[0] + r*np.cos(np.pi/3), origin[1] + r*np.sin(np.pi/3))
p5 = (origin[0] + r*np.cos(5*np.pi/3), origin[1] + r*np.sin(5*np.pi/3))
nodes = [p0,p1,p5]

G.add_node(i,pos=origin,center=True, phase_angle = np.random.randint(2000), boundary=[36,37,29,31,33,38])
i +=1
    
for node in nodes:
    G.add_node(i,pos=node,center=False)
    i += 1

###########################

    
# edges that are only passively elastic
G.add_path([1,2,3,4,5,6,1],beta=0,myosin=0, color = 'r')
G.add_path([1,17,16,8,9,2],beta=0,myosin=0, color = 'r')
G.add_path([6,11,12,18,17,1],beta=0,myosin=0,color = 'r')
G.add_path([18,19,14,15,16],beta=0,myosin=0, color = 'r')
G.add_path([15,21,22,23,8],beta=0,myosin=0,color='r')
G.add_path([12,26,27,25,19],beta=0,myosin=0,color='r')
G.add_path([25,34,33,31,14],beta=0,myosin=0,color='r')
G.add_path([31,29,30,21],beta=0,myosin=0,color='r')
G.add_path([33,38,36,37,29],beta=0,myosin=0,color='r')
    
# edges that have active force component
G.add_edges_from([(0,1),(0,2),(0,3),(0,4),(0,5),(0,6)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(7,16),(7,8),(7,9),(7,2),(7,1),(7,17)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(10,18),(10,17),(10,1),(10,6),(10,11),(10,12)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(13,14),(13,15),(13,16),(13,17),(13,18),(13,19)],beta=const.beta,myosin=const.myo0)

G.add_edges_from([(20,21),(20,22),(20,23),(20,8),(20,16),(20,15)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(24,25),(24,19),(24,18),(24,12),(24,26),(24,27)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(32,33),(32,31),(32,14),(32,19),(32,25),(32,34)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(28,29),(28,30),(28,21),(28,15),(28,14),(28,31)],beta=const.beta,myosin=const.myo0)
G.add_edges_from([(35,36),(35,37),(35,29),(35,31),(35,33),(35,38)],beta=const.beta,myosin=const.myo0)


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
dt = 0.1
tf = 6000
(t,Ac,Am,Bc,Bm,Rm,AB,AR) = dde_initializer(Ac,Am,Bc,Bm,Rm,AB,AR,tf,dt)


pic_num = 0 
for index in range(0,len(t)):
    H = G.copy() 
    ## Update myosin concentration on each spoke ##
    myo_list = []
    for center in G.nodes_iter(data=True):
        if center[1]['center']==True:
            # Calculate cell area
            outer = [G.node[element]['pos'] for element in center[1]['boundary']]
            cell_area = CellArea(outer)
            
            myosin_total = 0        # zero total myosin before continuing to next cell
            
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
                
                # Sum the total myosin in the current cell
                myosin_total += myosin_current

                #update myosin on this edge
                G[center[0]][neighbor]['myosin'] = dmyosin(myosin_current, geo_frac*Rm[index-center[1]['phase_angle']], length, dt)
            
#            print(center[0], myosin_total)
            myo_list.append(myosin_total)
    
#    print(myo_list)
    BioParamsWriter.writerow([t[index],myo_list[0],myo_list[1],myo_list[2],myo_list[3],myo_list[4],myo_list[5],myo_list[6],myo_list[7],myo_list[8]])


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

        if point not in [2,3,4,5,6,11,12,26,27,25,34,33,38,36,37,29,30,21,22,23,8,9]:
            G.node[point]['pos'] = d_pos(H.node[point]['pos'],total_force, dt)

#    BioParamsWriter.writerow([t[index], Rm[index], G[0][1]['myosin'], G[1][17]['myosin'], G[1][7]['myosin'], G[17][13]['myosin'], G[1][10]['myosin'],G[10][17]['myosin']])

    if index % 10 == 0:
        pic_num += 1
        print t[index]
        plt.clf()
        pos = nx.get_node_attributes(G,'pos')

        edges,colors = zip(*nx.get_edge_attributes(G,'color').items())
        nx.draw(G,pos, node_size = 1, edgelist=edges,edge_color=colors,width=1)
        
        plt.xlim(-20,70)
        plt.ylim(-25,25)
        plt.axis("on")
        plt.grid("on")
        plt.suptitle("t = %s"%t[index])

        plt.savefig('tmp%03d.png'%pic_num)


