""" Non-negative matrix tri-factorization (numpy)"""

import numpy as np
from scipy.linalg import svd
from math import sqrt

class SVD_EMBEDDINGS(object):
    def __init__(self,data_matrix):
        self.U, self.S, self.V = svd(data_matrix, full_matrices=True)

    "Compute SVD Based Embeddings"
    def embedding_space_svd(self,k):
        S = [sqrt(i) for i in self.S]
        diag_singular_values = np.diag(S)
        factorized_matrix = np.dot(self.U[:,0:k],diag_singular_values[0:k,0:k])
        return factorized_matrix


class SNMTF:
    """Compute Partial Orthonormal Non-negative Matrix Tri-Factorization (NMTF) Embeddings"""
    def __init__(self, max_iter=1000, verbose = 10):

        super(SNMTF,self).__init__()
        self.max_iter = max_iter
        self.verbose = verbose


    def Score(self, R1, P, U, G, norm_R1):
        GT = np.transpose(G)

        ErR1 = R1 - np.matmul(P, np.matmul(U, GT))
        norm_erR1 = np.linalg.norm(ErR1, ord='fro')

        rel_erR1 = norm_erR1 / norm_R1

        return norm_erR1, rel_erR1


    def Get_Clusters(self, M, nodes):
        n,k = np.shape(M)
        Clusters = [[] for i in range(k)]
        for i in range(n):
            idx = np.argmax(M[i])
            Clusters[idx].append(nodes[i])
        return Clusters



    def Solve_MUR(self, R1, k1, init="rand"):
        print "Starting NMTF"
        n,n=np.shape(R1)
        norm_R1 = np.linalg.norm(R1, ord='fro')

        if init == "rand":
            print " * Using random initialization"
            P = np.random.rand(n,k1) + 1e-14
            U = np.random.rand(k1,k1) + 1e-14
            G = np.random.rand(n,k1) + 1e-14
        elif init == 'SVD':
            # Initialize matrix factor P, U, and G using SVD decomposition

            print " * -- Eig decomposition on R1 to get PxUxG^T"

            J,K,L = svd(R1)

            P = np.zeros((n,k1)) + 1e-14
            G = np.zeros((n,k1)) + 1e-14
            U = np.zeros((k1,k1)) + 1e-14

            for i in range(n):
                for j in range(k1):
                    if J[i][j] > 0.:
                        P[i][j] = J[i][j]

            for i in range(k1):
                U[i][i] = abs(K[i])

            for i in range(n):
                for j in range(k1):
                    if L[i][j] > 0.:
                        G[i][j] = L[i][j]

        else :
            print "Unknown initializer: %s"%(init)
            exit(0)
        OBJ, REL2 = self.Score(R1,P, U, G, norm_R1)
        print " - Init:\t OBJ:%.4f\t REL2:%.4f"%(OBJ, REL2)


        #Begining M.U.R.
        for it in range(1,self.max_iter+1):

            R1T = np.transpose(R1)
            UT = np.transpose(U)

            PT = np.transpose(P)
            PT_P = np.matmul(PT,P)

            GT = np.transpose(G)
            GT_G = np.matmul(GT,G)

            #update rule for G

            #R1 matrix

            #update rule for P

            R1_G_UT = np.matmul(np.matmul(R1, G), UT)
            U_GT_G_UT = np.matmul(np.matmul(U,GT_G),UT)

            P_U_GT_G_UT = np.matmul(P,U_GT_G_UT)

            P_mult = np.sqrt( np.divide(R1_G_UT,P_U_GT_G_UT + 1e-4) )

            #update rule for U

            PT_R1_G = np.matmul(PT, np.matmul(R1, G))
            PT_P_U_GT_G = np.matmul(PT_P, np.matmul(U, GT_G))

            U_mult = np.sqrt( np.divide(PT_R1_G,PT_P_U_GT_G + 1e-4))

            #update rule for G

            R1T_P_U = np.matmul(np.matmul(R1T, P), U)
            G_GT_R1T_P_U = np.matmul(G, np.matmul(GT, R1T_P_U))

            G_mult = np.sqrt( np.divide(R1T_P_U, G_GT_R1T_P_U + 1e-14) )

            # Applying M.U.R.
            P = np.multiply(P, P_mult) + 1e-14
            U = np.multiply(U, U_mult) + 1e-14
            G = np.multiply(G, G_mult) + 1e-14


            if (it%self.verbose == 0) or (it==1):
                OBJ, REL2 = self.Score(R1, P, U, G, norm_R1)
                print " - It %i:\t OBJ:%.4f\t REL2:%.4f"%(it, OBJ, REL2)
        return P, U, G
