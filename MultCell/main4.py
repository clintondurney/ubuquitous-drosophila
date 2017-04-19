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
# Last Edit: 4/17/17
#
# To Do:    
#
#
###########

centers = [3,8,12,15,20,27,23,30,37,44,33,40,47,54,61,50,57,64,71,67,74,79,82,86,91]

# Set-up output file: myosin
myosinFile = open('myosin.csv','w')
myosinWriter = csv.writer(myosinFile,delimiter='\t')
header = ['cell' + str(node) for node in centers]
header = ['time'] + header
myosinWriter.writerow(header)

# Set-up output files: cell area
cellareaFile = open('cellarea.csv','w')
cellareaWriter = csv.writer(cellareaFile,delimiter='\t')
cellareaWriter.writerow(header)

# Set-up output files: cell area
regFile = open('reg.csv','w')
regWriter = csv.writer(regFile,delimiter='\t')
regWriter.writerow(header)

# Initialize tissue
G = tissue()
H = nx.Graph()

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
tf = 10000
(t,Ac,Am,Bc,Bm,Rm,AB,AR) = dde_initializer(Ac,Am,Bc,Bm,Rm,AB,AR,tf,dt)

pic_num = 0 
for index in range(0,len(t)):
    H = G.copy() 
    
    nodes = nx.get_node_attributes(G,'pos')
    
    ## Update myosin concentration on each spoke ##
    myo_hist = []
    area_hist = []
    reg_hist = []

    for n in centers:
        myosin_total = 0        # zero total myosin before continuing to next cell
        
        # Calculate cell area for the nth cell
        corners = [neighbor for neighbor in G.neighbors(n)]
        corn_sort = [(corners[0],0)]
        u = unit_vector(nodes[n],nodes[corners[0]])
        for i in range(1,len(corners)):
            v = unit_vector(nodes[n],nodes[corners[i]])
            dot = np.dot(u,v)
            det = np.linalg.det([u,v])
            angle = np.arctan2(det,dot)
            corn_sort.append((corners[i],angle))
        corn_sort = sorted(corn_sort, key=lambda tup: tup[1])
        corn2 = [nodes[entry[0]] for entry in corn_sort]
        
        cell_area = CellArea(corn2)
        
        for j in range(0,len(corn2)):
            # Calculate area of adjacent triangles of the spoke    
            inner = [corn2[np.mod(j,6)],corn2[np.mod(j+1,6)],nodes[n],corn2[np.mod(j-1,6)]]
             
            spoke_area = CellArea(inner)
            geo_frac = spoke_area/cell_area
            
            # Calculate necessary parameters for dm/dt
            length = distance.euclidean(nodes[n],corn2[j])
            myosin_current = G[n][corn_sort[j][0]]['myosin']
            
            #update myosin on this edge
            if (index - G.node[n]['time_lag']) >= 0:
                Reg = Rm[index-G.node[n]['time_lag']]
            else:
                Reg = 0
            G[n][corn_sort[j][0]]['myosin'] = dmyosin(myosin_current, geo_frac*Reg, length, dt)

            # Sum the total myosin in the current cell
            myosin_total += G[n][corn_sort[j][0]]['myosin']

        myo_hist.append(myosin_total)
        area_hist.append(cell_area)
        reg_hist.append(Reg)

    myosinWriter.writerow([t[index]] + myo_hist)
    cellareaWriter.writerow([t[index]] + area_hist)
    regWriter.writerow([t[index]] + reg_hist)

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

        if point not in [0,1,4,5,9,10,16,17,24,25,34,42,51,59,68,76,83,88,92,94,93,90,89,85,84,78,77,70,69,60,52,43,35,26,18,11,6,2]:
            G.node[point]['pos'] = d_pos(H.node[point]['pos'],total_force, dt)

    if index % 10 == 0:
        pic_num += 1
        print t[index]
        plt.clf()
        pos = nx.get_node_attributes(G,'pos')

        edges,colors = zip(*nx.get_edge_attributes(G,'color').items())
        nx.draw(G,pos, node_size = 1, edgelist=edges,edge_color=colors,width=1)
        
        plt.axis("on")
        plt.grid("on")
        plt.suptitle("t = %s"%t[index])

        plt.savefig('tmp%03d.png'%pic_num)


