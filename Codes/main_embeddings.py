import os
os.environ["OMP_NUM_THREADS"] = "4" # export OMP_NUM_THREADS=4
os.environ["OPENBLAS_NUM_THREADS"] = "4" # export OPENBLAS_NUM_THREADS=4
os.environ["MKL_NUM_THREADS"] = "4" # export MKL_NUM_THREADS=6
os.environ["VECLIB_MAXIMUM_THREADS"] = "4" # export VECLIB_MAXIMUM_THREADS=4
os.environ["NUMEXPR_NUM_THREADS"] = "4" # export NUMEXPR_NUM_THREADS=6

import networkx as nx
import numpy as np
import network_matrices as nm
import generate_embedding_space
import math
import time
import sys

print "\n\x1b[1;37;44m ############################################# \x1b[0m"
print "\x1b[1;37;44m #                                           # \x1b[0m"
print "\x1b[1;37;44m # Embeddings                                # \x1b[0m"
print "\x1b[1;37;44m #                                           # \x1b[0m"
print "\x1b[1;37;44m ############################################# \x1b[0m"

#Loading and preparing the matrix
type_of_matrix = sys.argv[1]

if type_of_matrix == 'ADJ':
    fname = "./Data/biogrid_homo_sapiens_ppi_edge_list"
    print "---- Loading PPI network"
    PPI, nodes, node2ind = nm.Load_Network(fname)
    PPI_mat = nm.Make_Adj(PPI, nodes, node2ind)
elif type_of_matrix == 'PPMI':
    fname = "./Data/biogrid_homo_sapiens_ppi_edge_list"
    print "---- Loading PPI network"
    PPI, nodes, node2ind = nm.Load_Network(fname)
    PPI_mat = nm.Make_Adj(PPI, nodes, node2ind)
    PPI_mat = nm.deep_walk_ppmi(PPMI_mat)
elif type_of_matrix == 'GDV':
    fname = './Data/gdvs_homo_sapiens_ppi.csv'
    print "---- Loading PPI network"
    PPI, nodes, node2ind = nm.Load_Network(fname)
    PPI_mat = nm.Make_Adj(PPI, nodes, node2ind)
elif type_of_matrix == 'GDV PPMI':
    fname = './Data/gdvs_homo_sapiens_ppi.csv'
    print "---- Loading PPI network"
    PPI, nodes, node2ind = nm.Load_Network(fname)
    PPI_mat = nm.Make_Adj(PPI, nodes, node2ind)
    PPI_mat = nm.deep_walk_ppmi(PPMI_mat)
else:
    raise NameError('Incorrect matrix type')

d = int(sys.argv[3])
features = [str(i) for i in range(d)]

decomposition_method = sys.argv[2]

if decomposition_method =='NMTF':
    Solver = generate_embedding_space.SNMTF(max_iter=500, verbose = 50)
    G, S = Solver.Solve_MUR(PPI_mat, d, init="SVD")
    nm.Save_Matrix_Factor(G, 'NMTF_'+type_of_matrix+'_'+d+'.csv', nodes, features)
else:
    SVD = generate_embedding_space.SVD_EMBEDDINGS(PPI_mat)
    G_SVD = embedding_space_svd(d)
    nm.Save_Matrix_Factor(G_SVD, 'SVD_'+type_of_matrix+'_'+d+'.csv', nodes, features)
