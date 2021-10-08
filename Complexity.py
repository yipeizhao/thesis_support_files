import numpy as np
import networkx as nx
from math import cos,pi,log
import utilities as ut
from scipy.linalg import eig
import math
#   All complexity measure have an optional parameter normalisation
#   If normalisation = True(set by default), the normalised value will be returned
#   Otherwise, the unnormalized form will be returned

#   Calculates OdC Complexity of a graph
def OdC(G,normalisation = True):
    if ut.empty_check(G) == True:
        return 0
    else:
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
            complexity = complexity/(ln(len(G.nodes)-1))
        return complexity

#   Calculates ln(x), if x=0, return 0
def ln(x):
    if x == 0:
        return 0
    else:
        return np.log(x)
    
#   Calculates Cr complexity of a graph
def Cr(G,normalisation = True):
    if ut.empty_check(G) == True:
        return 0
    else:
        if not nx.is_connected(G):
            return 0
        adj_matrix = nx.adjacency_matrix(G)
        adj_matrix = adj_matrix.todense()
        eigenvalues, eigenvectors = np.linalg.eig(adj_matrix)
        r = max(eigenvalues)
        r = r.real
        if normalisation == True:
            n = len(G.nodes)
            c_r_numerator = r-2 * cos(pi/(n+1))
            c_r_denominator = n-1-2*cos(pi/(n+1))
            c_r = c_r_numerator/c_r_denominator
            Cr_complexity = 4*c_r*(1-c_r)
            return Cr_complexity
        else:
            return r

#   Calculates Ce complexity of a graph
def Ce(G,normalisation = True):
    if not nx.is_connected(G):
        return 0
    if ut.empty_check(G):
        return -1
    else:
        E = 0
        nodes = list(G.nodes())
        n=len(nodes)
        for i in range(n):
            for j in range(i+1, n):
                E = E + 1/nx.shortest_path_length(G,nodes[i],nodes[j])
        E = 2* E /((n)*(n-1))
        if normalisation ==True:
            E_path = 0
            for i in range(1,n):
                E_path = E_path + (n-i)/i
            E_path = E_path *2/((n)*(n-1))
            Ce = 4*(E-E_path)/(1-E_path)*(1- (E-E_path)/(1-E_path))
            return Ce
        else:
            return E

#   Calculates C_{1e,st} complexity of a graph
def C1est(G,normalisation = True):
    n = len(G.nodes)
    mcu = n**1.68-10
    subgraphs = ut.subgraph_one_edge_deletion(G)
    st = []
    for item in subgraphs:
        L = nx.laplacian_matrix(item).toarray()
        L = np.delete(L,0,axis=0)
        L = np.delete(L,0,axis=1)
        _,logdet = np.linalg.slogdet(L)
        det = math.exp(logdet)
        det = int(det)
        det = round(det,10)
        if det not in st:
            st.append(det)
    N1est = len(st)
    if normalisation == False:
        return N1est
    else:
        return (N1est-1)/(mcu-1)

#   Calculates C_{1e.spec} complexity of a graph
def C1espec(G,normalisation =True):
    subgraphs = ut.subgraph_one_edge_deletion(G)
    spectra = [];spectra_s = []
    mcu = len(G.nodes())**1.68-10
    for i in range(len(subgraphs)):
        L = nx.laplacian_matrix(subgraphs[i]).todense()
        A = nx.adjacency_matrix(subgraphs[i]).todense()
        eig_values,_ = eig(L)
        eig_values = eig_values.real
        eig_values = sorted(eig_values)
        eig_values = [round(item,10) for item in eig_values]
        L = L + A + A
        eig_values_s,_ = eig(L)
        eig_values_s = eig_values_s.real
        eig_values_s = sorted(eig_values_s)
        eig_values_s = [round(item,10) for item in eig_values_s]
        if eig_values not in spectra and eig_values_s not in spectra_s:
            spectra.append(eig_values)
            spectra_s.append(eig_values_s)
    N1espec = len(spectra)
    if normalisation == False:
        return N1espec
    else:
        return (N1espec-1)/(mcu-1)

#   Calculates C_{2e,spec} complexity of a graph
def C2espec(G,normalisation=True):
    n= len(G.nodes)
    remove_edges = []
    for i in range(n):
        neighbours = list(G.neighbors(i))
        for item in neighbours:
            if item < i:
                remove_edges.append([i,item])
    products = []
    for i in range(len(remove_edges)):
        for j in range(i+1,len(remove_edges)):
            products.append([remove_edges[i],remove_edges[j]])
    subgraphs = []
    for item in products:
        temp_G = G.copy()
        temp_G.remove_edge(item[0][0],item[0][1])
        temp_G.remove_edge(item[1][0],item[1][1])
        subgraphs.append(temp_G)
    spectra = [];spectra_s = []
    mcu = len(G.nodes())**1.68-10
    for i in range(len(subgraphs)):
        L = nx.laplacian_matrix(subgraphs[i]).todense()
        A = nx.adjacency_matrix(subgraphs[i]).todense()
        eig_values,_ = eig(L)
        eig_values = eig_values.real
        eig_values = sorted(eig_values)
        eig_values = [round(item,10) for item in eig_values]
        L = L + A + A
        eig_values_s,_ = eig(L)
        eig_values_s = eig_values_s.real
        eig_values_s = sorted(eig_values_s)
        eig_values_s = [round(item,10) for item in eig_values_s]
        if eig_values not in spectra and eig_values_s not in spectra_s:
            spectra.append(eig_values)
            spectra_s.append(eig_values_s)
    N2espec = len(spectra)
    if normalisation == False:
        return N2espec
    else:
        C2espec = (N2espec - 1)/(math.comb(int(mcu),2))
        return C2espec

#   Calculates MA_{g} of a graph
#   Using function mutual_info and redundancy from utilities to calcualte I and R
def MAg(G,normalisation = True):
    n = len(G.nodes)
    R = ut.redundancy(G)
    I = ut.mutual_info(G)
    R_p = 2*(n-2)/(n-1)*log(2)
    R_c = 2*log(n-1)
    I_p = log(n-1)-((n-3)/(n-1))*log(2)
    I_c = log((n)/(n-1))
    MAr = 4*((R-R_p)/(R_c - R_p))*(1 - (R-R_p)/(R_c-R_p))
    MAi = 4*((I-I_c)/(I_p-I_c))*(1-(I-I_c)/(I_p-I_c))
    if normalisation == True:
        return MAr * MAi
    else:
        return R*I
        
#   Calcualtes MA_{RI} complexity of a graph
#   Using function mutual_info and redundancy from utilities to calcualte I and R
def MAri(G,normalisation=True):
    n = len(G.nodes)
    R = ut.redundancy(G)
    I = ut.mutual_info(G)
    R_p = 2*log(2)*(n-2)/(n-1)
    R_c = 2*log(n-1)
    I_p = log(n-1)-log(2)*(n-3)/(n-1)
    I_c = log(n/(n-1))
    m=len(G.edges)
    numerator_1 = (R-R_p)
    numerator_2 = (I-I_c)
    denominator = 0.25*(log(2*m)-R_p-I_c)**2
    if normalisation == True:
        return numerator_1*numerator_2/denominator
    else:
        return R*I