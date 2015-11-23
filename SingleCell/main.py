import numpy as np
from scipy.integrate import odeint
import pylab as plt
import networkx as nx
import globals as const
from dde_main import * 


############
#
# main.py
#
#
# Author: Clinton H. Durney
# Email: cdurney@math.ubc.ca
#
# Last Edit: 13/11/15
###########
    
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
G.add_node('medial',{'pos':origin, 'cell':1})

i = 0
for node in nodes:
    G.add_node(i,{'pos':node,'cell':1})
    G.add_edge('medial',i)
    i += 1
G.add_path([0,1,2,3,4,5,0])

# Initial plot of cell
pos = nx.get_node_attributes(G,'pos')
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

tf = 10 
(t,Ac,Am,Bc,Bm,Rm,AB,AR) = dde_solver(Ac,Am,Bc,Bm,Rm,AB,AR,tf)

tf = 1500 
(t,Ac,Am,Bc,Bm,Rm,AB,AR) = dde_tester(t,Ac,Am,Bc,Bm,Rm,AB,AR,tf)

tf = 2000
(t,Ac,Am,Bc,Bm,Rm,AB,AR) = dde_tester(t,Ac,Am,Bc,Bm,Rm,AB,AR,tf)

tf = 2500
(t,Ac,Am,Bc,Bm,Rm,AB,AR) = dde_tester(t,Ac,Am,Bc,Bm,Rm,AB,AR,tf)

tf = 6000
(t,Ac,Am,Bc,Bm,Rm,AB,AR) = dde_tester(t,Ac,Am,Bc,Bm,Rm,AB,AR,tf)

###
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
# G.node[0]['pos'] = (1,1)
# pos = nx.get_node_attributes(G,'pos')
# nx.draw(G,pos)
# plt.show()
###







