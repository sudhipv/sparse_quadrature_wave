

###### Trial



# from dolfin import *
import numpy as np
import time
import os as os
import scipy.io as sio
import numpy.random as npr


from mpi4py import MPI

## Initialize MPI for python
comm = MPI.COMM_WORLD

ip = comm.Get_rank()

if(ip ==0):
    t1 = time.time()

print('My rank is ',ip)


if(ip == 0):
    nSamples = 100
    NP = 4
    chunks = int(nSamples/NP)
    print('chunks is ', chunks)
    var_mu_sum = 0
    var_global = 0
else:
    chunks = None
    NP = None



var_local = 0


chunks = comm.bcast(chunks, root = 0)
NP = comm.bcast(NP, root = 0)


print('chunks is ', chunks)

# npr.seed(42)
xi = npr.randn(chunks)*10

print('xi values', xi, 'from',ip)


mean_local = sum(xi)/chunks

print('mean local is',mean_local)

for j in range(chunks):
    # print('xi[j] is', xi[j])
    var_local = var_local + np.power((xi[j] - mean_local),2)


var_local = var_local/(chunks-1)

print('var local is',var_local)

mean_sum = comm.allreduce(mean_local, op=MPI.SUM)

print('mean sum is', mean_sum)

var_global = comm.reduce(var_local, op=MPI.SUM, root=0)

if(ip ==0):
    print('var global is', var_global)

mu_exact = mean_sum/NP

var_mu = np.power((mean_local - mu_exact),2)

var_mu_sum = comm.reduce(var_mu, op=MPI.SUM, root=0)



if(ip == 0):
    sd_sum = var_global * (chunks-1)/(nSamples-1)

    v_mu = var_mu_sum * (chunks)/(nSamples-1)

    v_exact = sd_sum + v_mu

    print('mean of total samples is', mu_exact)
    print('variance of total samples is', v_exact)
    # sd_exact = np.sqrt(v_exact)






