import pdb
import csv
import numpy as np
import pylab as plt
import networkx as nx
import globals as const
from dde_main import * 
from scipy.spatial import distance

############
#
# main.py
#
#
# Author: Clinton H. Durney
# Email: cdurney@math.ubc.ca
#
# Last Edit: 03/05/16
#
# To Do:    
#           write loop instead of explicit calcs
#           Update loop to account for all spokes
#           Write other script to parse the output
#           Update spatial distribution of signal
#
#
###########

# initialize output files
LocationFile = open('Location.csv','w')
LocationWriter = csv.writer(LocationFile, delimiter='\t')
LocationWriter.writerow(["time","x_loc_0","y_loc_0","x_loc_1","y_loc_1"])

BioParamsFile = open('BioParams.csv','w')
BioParamsWriter = csv.writer(BioParamsFile,delimiter='\t')
BioParamsWriter.writerow(["time","myosin0","Force0","myosin1","Force1"])

num_cells = 1       # number of cells

G = nx.Graph()

# initial cell coordinates
origin = (0,0)
p0 = (1,0)
p1 = (np.cos(np.pi/3),np.sin(np.pi/3))
p2 = (np.cos(2*np.pi/3),np.sin(2*np.pi/3))
p3 = (np.cos(3*np.pi/3),np.sin(3*np.pi/3))
p4 = (np.cos(4*np.pi/3),np.sin(4*np.pi/3))
p5 = (np.cos(5*np.pi/3),np.sin(5*np.pi/3))
nodes = [p0,p1,p2,p3,p4,p5]

# initalize Graph G of nodes and edges
#G.add_node('medial',{'pos':origin, 'cell':1})
#
#i = 0
#for node in nodes:
#    G.add_node(i,{'pos':node,'cell':1})
#    G.add_edge('medial',i,{'name':1})
#    i += 1
#G.add_path([0,1,2,3,4,5,0])
#
## Initial plot of cell
#pos = nx.get_node_attributes(G,'pos')
#nx.draw(G,pos)
#plt.show()

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


# Initialize variables for the mechanical model
# Mar. 5 - As of now, I am appending numpy arrays so that they have a history to them.
#           If I am outputting to a csv file for interpreting later, I only need to keep 
#           the most recent values.

# location of nodes
p0_loc = [[p0[0],p0[1]]]
p1_loc = [[p1[0],p1[1]]]
p2_loc = [[p2[0],p2[1]]]
p3_loc = [[p3[0],p3[1]]]
p4_loc = [[p4[0],p4[1]]]
p5_loc = [[p5[0],p5[1]]]

# myosin conc. on spoke (variable)
myosin = np.array([[const.myo0,const.myo0]])
# initial length of spoke
l0 = np.array([[distance.euclidean(origin,p0),distance.euclidean(origin,p1)]])
# length of spoke
length = np.array([[distance.euclidean(origin,p0),distance.euclidean(origin,p1)]])
# force -- maybe should actually calculate initial force?
force = np.array([[0,0]])

# Time difference using discretization provided by the dde solver
dt = np.diff(t)

for index in range(0,len(dt)):
    # At each time step:
    if t[index] >= 0:
        
        # Update myosin concentration on each spoke
        temp = []
        for i in range(0,2):
            temp.append(dmyosin(myosin[-1][i],Rm[index], length[-1][i], dt[index]))
        myosin = np.append(myosin, [temp],axis=0)
        
        # Update force
        temp = []
        for i in range(0,2):
            temp.append(calc_force(length[-1][i], myosin[-1][i])) 
        force = np.append(force,[temp],axis=0)
    
        # Update Location
        # in order to loop over this, I need to put the node locations into one list.  Can't loop over the names...
        direction0 = normed_direction(origin,p0)
        direction1 = normed_direction(origin,p1)
        
        p0_loc.append(d_pos(p0_loc[-1][0],p0_loc[-1][1],force[-1][0],direction0,dt[index]))
        p1_loc.append(d_pos(p1_loc[-1][0],p1_loc[-1][1],force[-1][1],direction1,dt[index]))

        # update Length  
        length0 = distance.euclidean(p0_loc[-1],origin) 
        length1 = distance.euclidean(p1_loc[-1],origin) 
        length = np.append(length,[[length0,length1]],axis=0)

        # Write output file
        LocationWriter.writerow([t[index],p0_loc[-1][0],p0_loc[-1][1],p1_loc[-1][0],p1_loc[-1][1]])
        BioParamsWriter.writerow([t[index],myosin[-1][0], force[-1][0],myosin[-1][1],force[-1][1]])

        if p0_loc[-1][0] < 0:
            print index
            break


#############################################
#                                           #
# Plotting                                  #
#                                           #
#############################################

plt.figure(5)
plt.title("$Myosin$")
plt.plot(tt, myosin,'r')
plt.xlim([tt[0],tt[-1]])
plt.xlabel("Time (s)")
plt.ylabel("Myosin Concentration$")

plt.figure(6)
plt.title("Force")
plt.plot(tt,force)
plt.xlim([tt[0],tt[-1]])
plt.xlabel("Time (s)")
plt.ylabel("Force (N)")

plt.figure(7)
plt.title("Spoke Location")
plt.plot(tt,x_loc)
plt.xlim([tt[0],tt[-1]])
plt.xlabel("Time (s)")
plt.ylabel("Location")

# Create plots of biochemical parameters
plt.figure(1)
plt.title("$Reg_m$")
plt.plot(t, Rm,'r')
plt.xlim([t[0],t[-1]])
plt.xlabel("Time (s)")
plt.ylabel("$Reg_m$")

plt.figure(2)
plt.title("$aPKC_m$")
plt.plot(t,Am,'g')
plt.xlim([t[0],t[-1]])
plt.xlabel("Time (s)")
plt.ylabel("$aPKC_m$")

plt.figure(3)
plt.title("$Baz_m$")
plt.plot(t,Bm,'b')
plt.xlim([t[0],t[-1]])
plt.xlabel("Time (s)")
plt.ylabel("$Baz_m$")

plt.figure(4)
plt.title("Medial Reg, aPKC and Baz vs. Time")
plt.xlim([t[0],t[-1]])
plt.plot(t, Rm/np.amax(Rm),'r', label = "$Reg_m$")
plt.plot(t,Am/np.amax(Am), 'g', label = "$aPKC_m$")
plt.plot(t,Bm/np.amax(Bm), 'b', label = "$Baz_m$")
plt.legend()

plt.show()

# Example of manually changing the coordinates of point and replotting
#G.node[0]['pos'] = (0.75,0)
#G.node[1]['pos'] = (np.cos(np.pi/3)-.2,np.sin(np.pi/3)-.2) 
#pos = nx.get_node_attributes(G,'pos')
#nx.draw(G,pos)
#plt.show()
###


LocationFile.close()
BioParamsFile.close()


