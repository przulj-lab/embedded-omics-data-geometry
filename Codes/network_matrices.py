import numpy as np
from math import sqrt
import networkx as nx

def deep_walk_ppmi(adj_matrix, context_window=10):
    '''
    Input: Adj Matrix as numpy array
    '''
    degrees = np.sum(adj_matrix,axis = 0)
    volume_of_graph = sum(degrees)
    diag_degrees_inv = np.diag([1/float(i) for i in degrees])

    pmi = np.zeros([len(adj_matrix),len(adj_matrix)])

    transition_matrix = np.matmul(diag_degrees_inv,adj_matrix)
    pmi_temp = transition_matrix
    pmi += transition_matrix

    for r in range(1,context_window):
        print 'Iteration: ',r
        pmi_temp = np.matmul(pmi_temp,transition_matrix)
        pmi += pmi_temp

    pmi = pmi/float(context_window)
    pmi = volume_of_graph*np.matmul(pmi,diag_degrees_inv)

    pmi_matrix = np.log(pmi, out=np.zeros_like(pmi), where=(pmi!=0))

    for i in pmi_matrix:
        i[i<0]=0

    return pmi_matrix

def Make_Adj(net, net_nodes, net_n2i):
    nb_nodes = len(net_nodes)
    A = np.zeros((nb_nodes,nb_nodes))

    for e in net.edges():
        if e[0] in net_n2i and e[1] in net_n2i:
            n1=net_n2i[ e[0] ]
            n2=net_n2i[ e[1] ]
            A[n1][n2] = 1.
            A[n2][n1] = 1.

    return A

def Load_Network(fname):
	net = nx.read_edgelist(fname)
	nodes = [n for n in net.nodes()]
	nodeset = set(nodes)
	nb_nodes = len(nodes)
	nodes2ind = {}
	for n in range(nb_nodes):
		nodes2ind[nodes[n]]=n
	return net, nodes, nodes2ind

def Save_Matrix_Factor(M, fname, rownames, colnames):
    n,m = M.shape
    ofile = open(fname, 'w')
    #column header
    for i in range(m):
       ofile.write("\t%s"%(colnames[i]))
       ofile.write("\n")
    #rows
    for i in range(n):
        ofile.write("%s"%(rownames[i]))
    for j in range(m):
        ofile.write("\t%s"%(str(M[i][j])))
        ofile.write("\n")
    ofile.close()
