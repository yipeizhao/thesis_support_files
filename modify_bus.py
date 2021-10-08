#Remove all nodes with degree 2 in a graph, aka modified bus
def remove_degree_2(G1):
    G= G1.copy()
    flag = 1
    while flag == 1:
        node_len = len(G.nodes)
        remove_list = []
        for node in list(G.nodes):
            if G.degree(node) == 2:
                remove_list.append(node)
        
        for node in remove_list:
            neighbors = tuple(G.neighbors(node))
            if len(neighbors) == 2:
                G.remove_node(node)
                G.add_edge(neighbors[0],neighbors[1])
        if node_len == len(G.nodes):
            flag =0
    return G