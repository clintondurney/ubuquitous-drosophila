import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import itertools
from scipy.spatial import distance
import globals as const


def tissue():

    r = const.l_initial
    xx = np.arange(0,15*r,r)
    yy = np.arange(0,11*r*np.sqrt(3)*0.5,r*np.sqrt(3)*0.5)
    points = list(itertools.product(xx,yy))
    
    # shear the points
    m = np.sqrt(3)/3
    
    new_points = []
    for i in range(0,len(points)):
        x = points[i][0] + m*points[i][1]
        y = points[i][1]
        new_points.append((x,y))  
        
    points = new_points
    
    G = nx.Graph()
    for i in range(0,len(points)):
        G.add_node(i,pos=points[i])
    
    section_1 = [0,1,2,3,4,11,12,13,14,22,23,24,25,33,34,35,44,45,46,55,56,66,67,77,88]
    section_2 = [121,132,133,143,144,145,154,155,156,157]
    section_3 = [7,8,9,10,19,20,21,31,32,43]
    section_4 = [76,87,98,109,120,131,142,153,164,97,108,119,130,141,152,163,118,129,140,151,162,139,150,161,160]
    
    G.remove_nodes_from(section_1)
    G.remove_nodes_from(section_2)
    G.remove_nodes_from(section_3)
    G.remove_nodes_from(section_4)
    
    G = nx.convert_node_labels_to_integers(G,first_label=0)
    
    pos = nx.get_node_attributes(G,'pos')
    col = nx.get_edge_attributes(G,'color')

    return(G) 
