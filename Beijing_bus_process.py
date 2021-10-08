import pandas as pd
import numpy as np
import networkx as nx
import utilities as ut
df = pd.read_csv("bus/raw/beijing.csv")
stops = set(df["站点名称"])
lines = set(df["线路名称"])
lines_id = np.linspace(1,len(lines),len(lines))
lines_id = [int(item) for item in lines_id]
stops_id = np.linspace(1,len(stops),len(stops))
stops_id = [int(item) for item in stops_id]

line_dict = dict(zip(lines_id,lines))
stop_dict = dict(zip(stops_id,stops))
stops_id_conversion =[]
for stop in df["站点名称"]:
    stops_id_conversion.append(list(stop_dict.values()).index(stop))
stops_id_conversion = [stop+1 for stop in stops_id_conversion]
df["Stop_ID"] = stops_id_conversion



lines_list_seperate = []
for item in list(lines):
    temp_list = df.loc[df["线路名称"]==item]
    lines_list_seperate.append(temp_list)
    
G = nx.Graph()
G.add_nodes_from(stops_id)

for lines in lines_list_seperate:
    for i in range(len(lines)):
        for first, second in zip(lines["Stop_ID"], lines["Stop_ID"][1:]):
            G.add_edge(first,second)

#Remove self-loops
G.remove_edges_from(nx.selfloop_edges(G))

edges = list(G.edges())
source = [item[0] for item in edges]
target = [item[1] for item in edges]
data = {"source" : source,"target" : target}
processed_beijing_bus = pd.DataFrame(data = data)
processed_beijing_bus.to_csv("bus/processed/beijing.csv",index = False)

Gcc = ut.gcc(G)
edges = list(Gcc.edges())
source = [item[0] for item in edges]
target = [item[1] for item in edges]
data = {"source" : source,"target" : target}
processed_beijing_bus_simple = pd.DataFrame(data = data)
processed_beijing_bus_simple.to_csv("bus/simple/beijing_bus_simple.csv",index= False)