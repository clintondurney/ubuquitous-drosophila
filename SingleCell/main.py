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
# Last Edit: 03/11/16
#
# To Do:    
#           Update spatial distribution of signal
#           Add in edges - Force
#
###########

# initialize output files
LocationFile = open('Location.csv','w')
LocationWriter = csv.writer(LocationFile, delimiter='\t')
LocationWriter.writerow(["time","x_0","y_0","x_1","y_1","x_2","y_2","x_3","y_3","x_4","y_4","x_5","y_5","area"])

BioParamsFile = open('BioParams.csv','w')
BioParamsWriter = csv.writer(BioParamsFile,delimiter='\t')
BioParamsWriter.writerow(["time","myo0","myo1","myo2","myo3","myo4","myo5", "force0","force1","force2","force3","force4","force5"])

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

node_loc = [p0_loc,p1_loc,p2_loc,p3_loc,p4_loc,p5_loc]
# node_loc[node][history][x or y]

# myosin conc. on spoke (variable)
myosin = np.array([[const.myo0 for i in range(0,6)]])
# myosin[history][spoke]

# initial length of spoke
l0 = np.array([ [distance.euclidean(origin,node_loc[i][-1][:]) for i in range(0,6)] ] )
# l0[history][spoke]

# length of spoke
length = np.array( [ [distance.euclidean(origin,node_loc[i][-1][:]) for i in range(0,6) ] ] )
# length = [history][spoke]

# force -- maybe should actually calculate initial force?
force = np.array([[0 for i in range(0,6)]])
# force[history,spoke]

# Time difference using discretization provided by the dde solver
dt = np.diff(t)

for index in range(1,len(dt)):
    # At each time step:
    if t[index] >= 0:
        
        # Update myosin concentration on each spoke
        temp = [dmyosin(myosin[-1][i],Rm[index], length[-1][i], dt[index]) for i in range(0,6)]
        myosin = np.append(myosin, [temp],axis=0)
        
        # Update force
        temp = [calc_force(length[-1][i], myosin[-1][i]) for i in range(0,6)] 
        force = np.append(force,[temp],axis=0)
    
        # Update Location
        node_loc_temp = []
        for i in range(0,6):
            direction = normed_direction(origin,node_loc[i][-1][:])
            node_loc_temp.append(d_pos(node_loc[i][-1][0],node_loc[i][-1][1],force[-1][i],direction,dt[index]))
        
        for i in range(0,6):
            node_loc[i].append(node_loc_temp[i])

        # Calculate area for output file
        points = [(node_loc[i][-1][0],node_loc[i][-1][1]) for i in range(0,6)] 
        area = CellArea(points)

        # update Length  
        temp = [distance.euclidean(node_loc[i][-1][:],origin) for i in range(0,6)]
        length = np.append(length, [temp], axis=0)    

        # Write output file
        LocationWriter.writerow([t[index], 
            node_loc[0][-1][0], node_loc[0][-1][1], 
            node_loc[1][-1][0], node_loc[1][-1][1],
            node_loc[2][-1][0], node_loc[2][-1][1],
            node_loc[3][-1][0], node_loc[3][-1][1],
            node_loc[4][-1][0], node_loc[4][-1][1],
            node_loc[5][-1][0], node_loc[5][-1][1], area ] )

        BioParamsWriter.writerow([t[index], myosin[-1][0], myosin[-1][1],
            myosin[-1][2], myosin[-1][3], myosin[-1][4], myosin[-1][5], 
            force[-1][0], force[-1][1], force[-1][2], force[-1][3],
            force[-1][4], force[-1][5]])


#############################################
#                                           #
# Plotting                                  #
#                                           #
#############################################

#
## Create plots of biochemical parameters
#plt.figure(1)
#plt.title("$Reg_m$")
#plt.plot(t, Rm,'r')
#plt.xlim([t[0],t[-1]])
#plt.xlabel("Time (s)")
#plt.ylabel("$Reg_m$")
#
#plt.figure(2)
#plt.title("$aPKC_m$")
#plt.plot(t,Am,'g')
#plt.xlim([t[0],t[-1]])
#plt.xlabel("Time (s)")
#plt.ylabel("$aPKC_m$")
#
#plt.figure(3)
#plt.title("$Baz_m$")
#plt.plot(t,Bm,'b')
#plt.xlim([t[0],t[-1]])
#plt.xlabel("Time (s)")
#plt.ylabel("$Baz_m$")
#
#plt.figure(4)
#plt.title("Medial Reg, aPKC and Baz vs. Time")
#plt.xlim([t[0],t[-1]])
#plt.plot(t, Rm/np.amax(Rm),'r', label = "$Reg_m$")
#plt.plot(t,Am/np.amax(Am), 'g', label = "$aPKC_m$")
#plt.plot(t,Bm/np.amax(Bm), 'b', label = "$Baz_m$")
#plt.legend()
#
#plt.show()
#

LocationFile.close()
BioParamsFile.close()


