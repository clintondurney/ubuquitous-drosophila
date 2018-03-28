import pdb
import networkx as nx
import matplotlib
matplotlib.use('agg')
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
#
###########

# Initialize tissue
G, centers, boundary, AS_boundary = tissue()
H = nx.Graph()

area_blacklist = []
pic_num = 0 

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

# Set-up output files: reg level 
# regFile = open('reg.csv','w')
# regWriter = csv.writer(regFile,delimiter='\t')
# regWriter.writerow(header)

# Run dde solver for the Biochemical concentrations
dt = 0.1
tf = 14000 
(t,Rm) = dde_baz_cluster(const.k1,const.k2,const.k3,const.k4,const.k5,const.k6,const.k7,const.k8,const.qR,const.tau1,const.Ac0,const.Bc0)

for t_index in range(0,len(t)+1):
    H = G.copy() 
    pos = nx.get_node_attributes(G,'pos')
    
    myo_hist = []
    area_hist = []
#    reg_hist = []

    ## Update myosin concentration on each spoke ##
    for n in centers:
        myosin_total = 0        # zero total myosin before continuing to next cell
        if n not in area_blacklist:
            # Calculate cell area for the nth cell
            # Need to sort outer nodes in CW direction to calculate area
            corners_need_sorted = [neighbor for neighbor in G.neighbors(n)]
            sorted_corners,corn_sort_deg = sort_corners(corners_need_sorted,pos[n],pos)
            cell_area = CellArea(sorted_corners)
            
            if cell_area < 1.0:
                delete_cell = True
                for neighbor in G.neighbors(n):
                    if distance.euclidean(pos[n],pos[neighbor]) > 0.5:
                        delete_cell = False
                        break
                if delete_cell == True:          
                    # if cell area is less than 1 micron^2, then remove cell by adding to area_blacklist
                    print "deleting cell ", n  
                    area_blacklist.append(n)
                    # Contract the neighboring nodes of the removed cell
                    for node_removed in G.neighbors(n):
                        if node_removed in AS_boundary:
                            AS_boundary.remove(node_removed)
                            if n not in AS_boundary:
                                # only add in once
                                AS_boundary.append(n)
                            # re-sort AS boundary
                            co_ords, node_deg_sorted = sort_corners(AS_boundary,(0,0),pos)
                            AS_boundary = [node_deg_sorted[j][0] for j in range(0,len(node_deg_sorted))]
                    
                        G = nx.contracted_nodes(G,n,node_removed)
                    
                        if node_removed in area_blacklist:
                            # if in this list, it is a former center that was contracted
                            G.remove_node(node_removed)
                    G.node[n]['frozen'] = determine_freeze(G,pos,n,boundary)
                    cell_area = 0
            else:
		geo_frac_list = []
                for j in range(0,len(sorted_corners)):
                    # Calculate area of adjacent triangles of the spoke    
                    num_corners = len(sorted_corners)
                    
                    tri_1 = [pos[n], sorted_corners[np.mod(j-1,num_corners)],sorted_corners[np.mod(j,num_corners)]]
                    tri_2 = [pos[n], sorted_corners[np.mod(j,num_corners)],sorted_corners[np.mod(j+1,num_corners)]]
                    
                    spoke_area = CellArea(tri_1)/2.0 + CellArea(tri_2)/2.0
                   
                    geo_frac = spoke_area/cell_area
                    geo_frac_list.append(geo_frac)

                    # Calculate necessary parameters for dm/dt
                    length = distance.euclidean(pos[n],sorted_corners[j])
                    myosin_current = G[n][corn_sort_deg[j][0]]['myosin']
                
                    # Update myosin on this edge (need to *10 because t_index is in tenths of seconds)
                    if (t_index - G.node[n]['time_lag']*10) >= 0:
                        Reg = Rm[t_index-G.node[n]['time_lag']*10]
                    else:
                        Reg = 0
                    
                    G[n][corn_sort_deg[j][0]]['myosin'] = dmyosin(myosin_current, geo_frac*Reg, length, dt)

                    # Sum the total myosin in the current cell
                    myosin_total += G[n][corn_sort_deg[j][0]]['myosin']
                
#                if myosin_total >= 40000:
#                    myosin_total = 40000
#                    for j in range(0,len(sorted_corners)):
#                        G[n][corn_sort_deg[j][0]]['myosin'] = geo_frac_list[j]*myosin_total
                
        else:
            cell_area = 0
        
        # Update list for CSV file writing
        myo_hist.append(myosin_total)
        area_hist.append(cell_area)
#       reg_hist.append(Reg)

    # Write to the CSV file
    myosinWriter.writerow([t[t_index]] + myo_hist)
    cellareaWriter.writerow([t[t_index]] + area_hist)
#    regWriter.writerow([t[t_index]] + reg_hist)


    ## Actin Cable Formation ##
    # Formation follows a Logisitic Curve with params AC_max, k (steepness) and x0 - center of log. curve
    for j in range(0,len(AS_boundary)):
        G[AS_boundary[j-1]][AS_boundary[j]]['myosin'] = const.AC_max/(1+np.exp((-const.k)*(t[t_index]-const.x0)))


    ## Update force ##
    # iterate over all nodes in graph
    for point in G.nodes_iter():
    	# iterate over all neighbors of node
    	update_location = True
        total_force = [0,0]
    	for neighbor in G.neighbors(point):
            if neighbor not in boundary:
                # if neighbor is not in boundary, then have passive and active forces
                # calculate the unit vector from node to neighbor
                dir_vector = unit_vector(H.node[point]['pos'],H.node[neighbor]['pos'])
                # calculate magnitude of force
                length = distance.euclidean(H.node[point]['pos'],H.node[neighbor]['pos'])
                mag_force = calc_force(length,G[point][neighbor]['myosin'])
            #    total_force = np.sum([total_force,mag_force*np.array(dir_vector)],axis=0)
            else:
                # if neighbor is in boundary, then it is an epidermis spoke
                # calculate the unit vector from node to neighbor
                dir_vector = unit_vector(H.node[point]['pos'],H.node[neighbor]['pos'])
                length = distance.euclidean(H.node[point]['pos'],H.node[neighbor]['pos'])
                # calculate magnitude of force
                # desired constant force -- F=mu*l 
                const_force_length = const.epi_tension/const.mu
                # mag_force = const. force + elastic force component
                mag_force = calc_force(const_force_length,0) + calc_force(length,0) 
            total_force = np.sum([total_force,mag_force*np.array(dir_vector)],axis=0)
        
        if G.node[point]['frozen'] == True:
            update_location = False

        # Update Node locations of those not fixed (the epidermis boundary)
        if update_location == True:
            G.node[point]['pos'] = d_pos(H.node[point]['pos'],total_force, dt)

    # Output a picture every 1 seconds
    if t_index % 10 == 0:
        pic_num += 1
        print t[t_index]
        plt.clf()
        pos = nx.get_node_attributes(G,'pos')

        edges,colors = zip(*nx.get_edge_attributes(G,'color').items())
        nx.draw(G,pos, node_size = 0, edgelist=edges,edge_color=colors,width=2,node_color='black')
        
        plt.axis("on")
        plt.grid("on")
        plt.axis("equal")
        plt.suptitle("t = %s"%t[t_index])

        plt.savefig('tmp%03d.png'%pic_num)


