


import numpy as np
import time
import os as os
import scipy.io as sio

from scipy.io import loadmat

from numpy.polynomial.hermite import Hermite



# def OnedPC(dim_L, xi):


#     PC_1d = np.zeros([11,dim_L])


#     PC_1d[0,:] = 1
#     PC_1d[1,:] = xi
#     PC_1d[2,:] = pow(xi,2) - 1
#     PC_1d[3,:] = pow(xi,3) - (3*xi)
#     PC_1d[4,:] = pow(xi,4) - (6*pow(xi,2)) + 3
#     PC_1d[5,:] = pow(xi,5) - (10*pow(xi,3)) + 15
#     PC_1d[6,:] = pow(xi,6) - (15*pow(xi,4)) + 45*pow(xi,2) -15
#     PC_1d[7,:] = pow(xi,7) - (21*pow(xi,5)) + 105*pow(xi,3) -105
#     PC_1d[8,:] = pow(xi,8) - (28*pow(xi,6)) + 210*pow(xi,4) - 420*pow(xi,2) + 105
#     PC_1d[9,:] = pow(xi,9) - (36*pow(xi,7)) + 378*pow(xi,5) - 1260*pow(xi,3) + 945*xi
#     PC_1d[10,:] = pow(xi,10) - (45*pow(xi,8)) + 630*pow(xi,6) - 3150*pow(xi,4) + 4725*pow(xi,2) - 945


#     return PC_1d



def recursivepc1D(index, xi):


    dim_L = len(xi)

    print(dim_L)

    PC_1d = np.zeros([index+1,dim_L])

    # print(PC_1d)

    if (index >= 2):

        PC_1d[0,:] = 1
        # print('index 0', PC_1d)

        PC_1d[1,:] = xi
        # print('index 1', PC_1d)

        for j in range(2, index+1):
            # print('j is', j)
            PC_1d[j,:] = np.multiply(xi, PC_1d[j-1,:]) - ((j-1) * PC_1d[j-2,:])
            # print('index ',j, PC_1d[j,:])

    elif (index >= 1 ):
        PC_1d[0,:] = 1
        PC_1d[1,:] = xi
        # print('index 1', PC_1d)

    else:
        PC_1d[0,:] = 1
        # print('index 0', PC_1d)

    return PC_1d



def multiPC(dim_L,PCE_u, mIndex, xi):

    PCmulti_quad = np.zeros([PCE_u])
    locindex = np.zeros([0,dim_L])

    tol = 1e-14

    for i in range(PCE_u):


        locindex = mIndex[i,:]

        print('locindex is', locindex)

        print('xi is', xi)

        index = max(locindex)

        print('max loc index is', index)

        PC1D = recursivepc1D(index, xi)

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


# PC1D = OnedPC(dim_L, xi)

PC_quad = multiPC(dim_L,PCE_u, mIndex, xi)



print(PC_quad)














# k = np.zeros([2,2])

# k[0,0] = 8
# k[0,1] = 6

# k[1,0] = 2
# k[1,1] = 4

# print(k)

# a = 5

# print(a)

# print('k *a is ',k*a)
