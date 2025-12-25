


import numpy as np
import time
import os as os
import scipy.io as sio
import math

from scipy.io import loadmat

from numpy.polynomial.hermite import Hermite



def nfactorial(nf):
    if nf == 0:
        return 1
    else:
        return nf * nfactorial(nf-1)

def OnedPC(dim_L, xi):


    PC_1d = np.zeros([11,dim_L])


    PC_1d[0,:] = 1
    PC_1d[1,:] = xi
    PC_1d[2,:] = pow(xi,2) - 1
    PC_1d[3,:] = pow(xi,3) - (3*xi)
    PC_1d[4,:] = pow(xi,4) - (6*pow(xi,2)) + 3
    PC_1d[5,:] = pow(xi,5) - (10*pow(xi,3)) + 15
    PC_1d[6,:] = pow(xi,6) - (15*pow(xi,4)) + 45*pow(xi,2) -15
    PC_1d[7,:] = pow(xi,7) - (21*pow(xi,5)) + 105*pow(xi,3) -105
    PC_1d[8,:] = pow(xi,8) - (28*pow(xi,6)) + 210*pow(xi,4) - 420*pow(xi,2) + 105
    PC_1d[9,:] = pow(xi,9) - (36*pow(xi,7)) + 378*pow(xi,5) - 1260*pow(xi,3) + 945*xi
    PC_1d[10,:] = pow(xi,10) - (45*pow(xi,8)) + 630*pow(xi,6) - 3150*pow(xi,4) + 4725*pow(xi,2) - 945


    return PC_1d



# def recursivepc1D(index, xi):


#     dim_L = len(xi)

#     print(dim_L)

#     PC_1d = np.zeros([index+1,dim_L])

#     # print(PC_1d)

#     if (index >= 2):

#         PC_1d[0,:] = 1
#         # print('index 0', PC_1d)

#         PC_1d[1,:] = xi
#         # print('index 1', PC_1d)

#         for j in range(2, index+1):
#             # print('j is', j)
#             PC_1d[j,:] = np.multiply(xi, PC_1d[j-1,:]) - ((j-1) * PC_1d[j-2,:])
#             # print('index ',j, PC_1d[j,:])

#     elif (index >= 1 ):
#         PC_1d[0,:] = 1
#         PC_1d[1,:] = xi
#         # print('index 1', PC_1d)

#     else:
#         PC_1d[0,:] = 1
#         # print('index 0', PC_1d)

#     return PC_1d
def normpc(dim_L,PCE_u, mIndex):

    normPC = np.zeros([PCE_u])

    tol = 1e-14

    for i in range(PCE_u):


        locindex = mIndex[i,:]

        PCmulti = 1
        normPsi = 1

        for d in range(dim_L):

            print('d is', d)

            xi_index = locindex[d]

            print('xi_index is',xi_index)

            normPsi = normPsi * nfactorial(xi_index)
            # print('PCmulti is',PCmulti)


        normPC[i] = normPsi

    return normPC


def multiPC(dim_L,PCE_u, mIndex, PC1D, xi):

    PCmulti_quad = np.zeros([PCE_u])
    locindex = np.zeros([0,dim_L])

    tol = 1e-14

    for i in range(PCE_u):


        locindex = mIndex[i,:]

        print('locindex is', locindex)

        print('xi is', xi)

        index = max(locindex)

        print('max loc index is', index)

        # PC1D = recursivepc1D(index, xi)

        print('PC 1D is', PC1D)

        PCmulti = 1

        for d in range(dim_L):

            print('d is', d)

            xi_index = locindex[d]

            print('xi_index is',xi_index)

            PCmulti = PCmulti * PC1D[xi_index,d]

            # print('PCmulti is',PCmulti)


        if(abs(PCmulti) < tol):
            PCmulti = 0

        PCmulti_quad[i] = PCmulti
        print('PCmulti 1 chaos term', PCmulti)

    return PCmulti_quad



dim_L = 3

PCE_u = 20

quadpath = "./quadData/d3l3/qdpts.dat"
weightpath = "./quadData/d3l3/wghts.dat"

qp = np.genfromtxt(quadpath)

norm_psi = loadmat('./klePceData/norm_squared030003.mat')

norm = norm_psi['norm_squared']

print('norm psi ', np.sqrt(norm[0,10]))


mIndexPath = "./klePceData/multiIndex.dat"
mIndex = np.genfromtxt(mIndexPath).astype(int)

xi = qp[23]



# if(ip ==0):
NP = 64

NQd = len(qp)

print('NP is', NP)

print('nqd is',NQd)

serpart = int(NQd/NP)

print('serpart is', serpart)

serpart = math.ceil(NQd/NP)

print('ceiled serpart is', serpart)

part1 = NP

part2 = NQd%NP

print('part2 is',part2)

chunks = len(qp)

print('chunks is', chunks)

# PC1D = OnedPC(dim_L, xi)

# PC_quad = multiPC(dim_L,PCE_u, mIndex,PC1D, xi)

# normPC = normpc(dim_L,PCE_u, mIndex)

# print(PC_quad)

# print(normPC)



############# EXTRA IMPORTANT STUFFS ########

        # if self.ndim == 2:


        #     #### Multiplying chaos terms with the Gaussian Sample for MCS
        #     xi_1 = self.xi[0]
        #     xi_2 = self.xi[1]
        #     # print('xi1 and xi2 is', xi_1,xi_2)


        #     lterms[0] = 1.0
        #     lterms[1] = xi_1
        #     lterms[2] = xi_2
        #     lterms[3] = pow(xi_1,2) - 1;
        #     lterms[4] = xi_1*xi_2;
        #     lterms[5] = pow(xi_2,2) - 1;
        #     lterms[6] = pow(xi_1,3) - 3*xi_1;
        #     lterms[7] = xi_1*xi_1*xi_2 - xi_2;
        #     lterms[8] = xi_2*xi_2*xi_1 - xi_1;
        #     lterms[9] = pow(xi_2,3) - 3*xi_2;


        # if self.ndim == 3:


        #     #### Multiplying chaos terms with the Gaussian Sample for MCS
        #     xi_1 = self.xi[0]
        #     xi_2 = self.xi[1]
        #     xi_3 = self.xi[2]
        #     # print('xi1 and xi2 is', xi_1,xi_2)

        #     lterms[0] = 1.0
        #     lterms[1] = xi_1;
        #     lterms[2] = xi_2;
        #     lterms[3] = xi_3;
        #     lterms[4] = pow(xi_1,2) - 1;
        #     lterms[5] = xi_1*xi_2;
        #     lterms[6] = xi_1*xi_3;
        #     lterms[7] = pow(xi_2,2) - 1;
        #     lterms[8] = xi_2*xi_3;
        #     lterms[9] = pow(xi_3,2) - 1;
        #     lterms[10] = pow(xi_1,3) - 3*xi_1;
        #     lterms[11] = pow(xi_1,2)*xi_2 - xi_2;
        #     lterms[12] = pow(xi_1,2)*xi_3 - xi_3;
        #     lterms[13] = pow(xi_2,2)*xi_1 - xi_1;
        #     lterms[14] = xi_1*xi_2*xi_3;
        #     lterms[15] = pow(xi_3,2)*xi_1 - xi_1;
        #     lterms[16] = pow(xi_2,3) - 3*xi_2;
        #     lterms[17] = pow(xi_2,2)*xi_3 - xi_3;
        #     lterms[18] = pow(xi_3,2)*xi_2 - xi_2;
        #     lterms[19] = pow(xi_3,3) - 3*xi_3;











# k = np.zeros([2,2])

# k[0,0] = 8
# k[0,1] = 6

# k[1,0] = 2
# k[1,1] = 4

# print(k)

# a = 5

# print(a)

# print('k *a is ',k*a)
