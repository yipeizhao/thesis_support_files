import numpy as np
import networkx as nx
from random import randint
from math import log
from scipy.linalg import eig
def C1espec(G,normalisation =True):
    subgraphs = subgraph_one_edge_deletion(G)
    spectra = []
    #Upper bound
    mcu = len(G.nodes())**1.68-10
    for i in range(len(subgraphs)):
        #L = laplacian matrix
        # L+2A to convert negative 1s to positive 1
        L = nx.laplacian_matrix(subgraphs[i]).todense()
        A = nx.adjacency_matrix(subgraphs[i]).todense()
        L = L + A + A
        #find the max eig value of L
        eig_values,_ = eig(L)
        spectra.append(max(eig_values.real))
    rounded_spectra = []
    #Round all the spectras, taking only the first 7 sf.
    for item in spectra:
        rounded_spectra.append(round(item,7))
    #Find the number of unique spectra values
    # #of unique spectra values = number of non-isomorphic subgraphs
    N1espec = len(set(rounded_spectra))
    if normalisation == False:
        return N1espec
    else:
        return (N1espec-1)/(mcu-1)

#Find all subgraph caused by deleting one edge
def subgraph_one_edge_deletion(G):
    subgraphs = []
    for edge in list(G.edges):
        temp_graph = nx.Graph(G)
        temp_graph.remove_edge(edge[0],edge[1])
        subgraphs.append(temp_graph)
    return subgraphs


def mutual_info(G):
    edges = list(G.edges)
    m = len(edges)
    I = 0
    for item in edges:
        d0 = len(list(G.neighbors(item[0])))
        d1 = len(list(G.neighbors(item[1]))) 
        I = I + log(2*m/(d0*d1))
    return I/m

def redundancy(G):
    edges = list(G.edges)
    m = len(edges)
    R = 0
    for item in edges:
        d0 = len(list(G.neighbors(item[0])))
        d1 = len(list(G.neighbors(item[1]))) 
        R = R + log((d0*d1))
    return R/m
def MAri(G,normalisation = True):
    n = len(G.nodes)
    R = redundancy(G)
    I = mutual_info(G)
    R_p = 2*log(2)*(n-2)/(n-1)
    R_c = 2*log(n-1)
    I_p = log(n-1)-log(2)*(n-3)/(n-1)
    I_c = log(n/(n-1))
    if normalisation == True:
        return 4 *((R-R_p)/(R_c - R_p))*((I-I_c)/(I_p-I_c))
    else:
        return R*I
def OdC(G,normalisation = True):
    #Create a degree correlation matrix, using the max degree 
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    max_degree = max(degree_sequence)
    degree_correlation = np.zeros((max_degree,max_degree))
    #Building the correlation matrix
    for node in list(G.nodes):
        #An array to store all the neighbors degrees
        neighbors_degree = []
        #Getting the degree of the current node
        node_degree = G.degree(node)
        #Stating all neighbors and finding their degrees
        neighbors = list(G.neighbors(node))
        neighbors_degrees_tuple=G.degree(neighbors)
        for item in neighbors_degrees_tuple:
            neighbors_degree.append(item[1])
        #For every occurence, adding one to the matrix
        for item in neighbors_degree:
            if node_degree<=item:
                degree_correlation[node_degree-1,item-1] +=1
    #Calculating a_k
    a_k=[]
    for i in range(max_degree):
        a_k.append(sum(degree_correlation[i]))
    A = sum(a_k)
    if A !=0:
        for i in range(len(a_k)):
            a_k[i]=a_k[i]/A
    #Calculating the complexity
    complexity = 0
    for item in a_k:
        complexity -= item*ln(item)
    #Normalisation
    if normalisation == True:
        complexity = complexity/(ln(len(G.degree)-1))
    return complexity

def ln(x):
    if x == 0:
        return 0
    else:
        return np.log(x)

def power_law_series(n,gamma):
    n = int(n)
    sequence = nx.random_powerlaw_tree_sequence(n, gamma,tries=100000)
    return sequence

def construct_network(edge_list):
    G=nx.Graph()
    for item in edge_list:
        G.add_edge(item[0],item[1])
    return G
