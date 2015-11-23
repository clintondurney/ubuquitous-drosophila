import networkx as nx
import pylab as plt
import numpy as np

##########
#
# cellnetwork.py
# keeping this separate for now (as a document of something that works), 
# but it has now been implemented into main.py
#
# Author: Clinton H. Durney
# Email: cdurney@math.ubc.ca
#
# Last Edit: 13/11/15
#########   

def plot(origin,p0,p1,p2,p3,p4,p5):
    nodes = [p0,p1,p2,p3,p4,p5]

    G = nx.Graph()

    G.add_node('medial',{'pos':origin, 'cell':1})
    i = 0
    for node in nodes:
        G.add_node(i,{'pos':node,'cell':1})
        G.add_edge('medial',i)
        i += 1

    G.add_path([0,1,2,3,4,5,0])

    pos = nx.get_node_attributes(G,'pos')

    nx.draw(G,pos)
    plt.show()



