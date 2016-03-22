# import csv
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


# Read into Pandas Data Frame
df_bio = pd.read_csv('BioParams.csv',delimiter='\t')
df_location = pd.read_csv('Location.csv',delimiter='\t')

# get the x and y positions into a list
time = df_bio['time']
x0 = df_location.x_0
y0 = df_location.y_0
x1 = df_location.x_1
y1 = df_location.y_1
x2 = df_location.x_2
y2 = df_location.y_2
x3 = df_location.x_3
y3 = df_location.y_3
x4 = df_location.x_4
y4 = df_location.y_4
x5 = df_location.x_5
y5 = df_location.y_5

p0 = (x0[0],y0[0])
p1 = (x1[0],y1[0])
p2 = (x2[0],y2[0])
p3 = (x3[0],y3[0])
p4 = (x4[0],y4[0])
p5 = (x5[0],y5[0])

nodes = [p0,p1,p2,p3,p4,p5]

G = nx.Graph()

i = 0
for node in nodes:
    G.add_node(i,{'pos':node,'cell':1})
#    G.add_edge('medial',i,{'name':1})
    i += 1
G.add_path([0,1,2,3,4,5,0])

fig = plt.figure()

pos = nx.get_node_attributes(G,'pos')

nx.draw(G,pos)
plt.title("Time")
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.axis("on")
plt.grid("on")

for i in range(1,len(time),1000):
	G.node[0]['pos'] = (x0[i],y0[i])
	G.node[1]['pos'] = (x1[i],y1[i])
	G.node[2]['pos'] = (x2[i],y2[i])
	G.node[3]['pos'] = (x3[i],y3[i])
	G.node[4]['pos'] = (x4[i],y4[i])
	G.node[5]['pos'] = (x5[i],y5[i])
	
	pos = nx.get_node_attributes(G,'pos')
	nx.draw(G,pos)
	plt.title('Time')
	plt.xlim(-2,2)
	plt.ylim(-2,2)
	plt.axis("on")
	plt.grid("on")
	plt.pause(0.5)

plt.show()
