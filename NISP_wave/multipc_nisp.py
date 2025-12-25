

####### Code to generate pdf from the NISP coefficients

import numpy as np
import time
import os as os
import scipy.io as sio
import numpy.random as npr
import pyvista as pv
import glob


def OnedPC(dim_L, xi):


    PC_1d = np.zeros([11,dim_L])

    # print('dim_L', dim_L)
    # print('xi inside onepc', xi)

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

#########################################################################


def normpc(dim_L,PCE_u, mIndex):

    normPC = np.zeros([PCE_u])

    tol = 1e-14

    for i in range(PCE_u):


        locindex = mIndex[i,:]

        normPsi = 1

        for d in range(dim_L):

            # print('d is', d)

            xi_index = locindex[d]

            # print('xi_index is',xi_index)

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


        PCmulti = 1

        for d in range(dim_L):

            # print('d is', d)

            xi_index = locindex[d]

            # print('xi_index is',xi_index)

            PCmulti = PCmulti * PC1D[xi_index,d]

            # print('PCmulti is',PCmulti)


        if(abs(PCmulti) < tol):
            PCmulti = 0

        PCmulti_quad[i] = PCmulti
        # print('PCmulti 1 chaos term', PCmulti)

    return PCmulti_quad



########### MAIN CODE ################

pcedatapath = "./klePceData/pcedata.dat"
pcedata = np.genfromtxt(pcedatapath)
order_u = pcedata[0].astype(int)     ## order of output PCE (input==2)
dim_L   = pcedata[1].astype(int)     ## number of KLE modes (RVs)
PCE_u   = pcedata[2].astype(int)     ## number of PCE outputs
PCE_A   = pcedata[3].astype(int)    ## number of PCE inputs

mIndexPath = "./klePceData/multiIndex.dat"
mIndex = np.genfromtxt(mIndexPath).astype(int)


print('Random Dimension', dim_L)
print('order of PC output', order_u)
print('number of PC output', PCE_u)



mainpath = '/Users/sudhipv/documents/python_ni/NISP_wave/results/'

path = '3rv_d3l3_order4_pdf/'

pc_coeff = np.zeros([9,PCE_u])
# print(pc_coeff)
matfile = mainpath + path + 'u_pdf.mat'


lfile = sio.loadmat(matfile)

######### Generating Multi dimensional PCs for each random variable to get PDF

pc_coeff = lfile['u_pdf']

samples = 300000

u_pdf = np.zeros([9,samples])


for i in range(samples):


    xi = npr.randn(dim_L)

    # print('xi is',xi)

    PC1D = OnedPC(dim_L, xi)

    lterms = np.zeros([PCE_u])

    lterms = multiPC(dim_L,PCE_u, mIndex, PC1D, xi)



    for k in range(9):

        for n_pc in range(PCE_u):

            u_pdf[k,i] = u_pdf[k,i] + pc_coeff[k,n_pc]*lterms[n_pc]


spath = mainpath + path + 'nisp_pdf_1'

sio.savemat(spath, {'nisp_pdf_1':u_pdf})


exit()
