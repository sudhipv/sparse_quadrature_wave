


### Acoustic Wave Propagation Problem - Serial ######

##### Monte Carlo Simulation using Fortran Decomposed Mesh and Fenics Assembly

#### Copyright (C) Sudhi P V #####

### June -30 - 2020



######## Variance Computation in parallel is tricky....Here its assumed that the sample sizes are same for all processes..if not the formulae changes..######

##################

from dolfin import *
import numpy as np
import time
import os as os
import scipy.io as sio
import numpy.random as npr


from mpi4py import MPI

## Initialize MPI for python
comm = MPI.COMM_WORLD

ip = comm.Get_rank()

## print('My rank is ',ip)

# if ip==0:
#     print("==========================================================")
#     print("Running FEniCS in Parallel for 2D stochastic Acoustic Assembly...")

##### VERY IMPORTANT - TO AVOID FENICS FROM RE ORDERING THE NODES
parameters['reorder_dofs_serial'] = False

## Global deterministic assembly for each KLE mode (here it's only mean term)
def detAssembly(a,m,l):

    # if (ip == 0):
    #     print("Inside Deterministic Assembly")
    ## Dummy problem to define LAYOUT of global problem (the easy way)
    A_g = assemble(Constant(0.)* c * inner(nabla_grad(u1),nabla_grad(w))*dx)
    M_g = assemble(Constant(0.)*u1*w*dx)
    L_g = assemble(Constant(0.)*f*w*dx)


    ## Get dofmap to construct cell-to-dof connectivity
    dofmap = V.dofmap()

    ## Perform assembling
    for cell in cells(mesh):
        dof_idx = dofmap.cell_dofs(cell.index())

        ## Assemble local rhs and lhs
        a_local  = assemble_local(a, cell)
        m_local  = assemble_local(m, cell)
        l_local  = assemble_local(l, cell)

        ## Assemble into global system
        A_g.add_local(a_local,dof_idx, dof_idx)
        M_g.add_local(m_local,dof_idx, dof_idx)
        L_g.add_local(l_local,dof_idx)

    ## Finalize assembling
    A_g.apply("add"), M_g.apply("add"), L_g.apply("add")

#################################################
    # if (k == 0):
    #     print("Setting boundary conditions")

    def boundary_S(x, on_boundary):
        tol = 1E-14
        return on_boundary


    # Mark boundary subdomians (2D domain (0,1) x (0,1)
    left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
    right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
    top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0)
    bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)

    bc_L = DirichletBC(V, Constant(0), left)
    bc_R = DirichletBC(V, Constant(0), right)
    bc_T = DirichletBC(V, Constant(0), top)
    bc_B = DirichletBC(V, Constant(0), bottom)

    ## Apply boundary conditions
    bcs = [bc_L,bc_R,bc_T,bc_B]


    for bc in bcs:
        bc.apply(A_g)
        bc.apply(M_g)
        bc.apply(L_g)

    A = A_g
    M = M_g
    b = L_g

    return A,M,b





def detAssemblyA(a):

    # if (ip == 0):
    #     print("Inside Deterministic Assembly")
    ## Dummy problem to define LAYOUT of global problem (the easy way)
    A_g = assemble(Constant(0.)* c * inner(nabla_grad(u1),nabla_grad(w))*dx)


    ## Get dofmap to construct cell-to-dof connectivity
    dofmap = V.dofmap()

    ## Perform assembling
    for cell in cells(mesh):
        dof_idx = dofmap.cell_dofs(cell.index())

        ## Assemble local rhs and lhs
        a_local  = assemble_local(a, cell)


        ## Assemble into global system
        A_g.add_local(a_local,dof_idx, dof_idx)


    ## Finalize assembling
    A_g.apply("add")

#################################################
    # if (k == 0):
    #     print("Setting boundary conditions")

    def boundary_S(x, on_boundary):
        tol = 1E-14
        return on_boundary


    # Mark boundary subdomians (2D domain (0,1) x (0,1)
    left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = 0.0)
    right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0)
    top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0)
    bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)

    bc_L = DirichletBC(V, Constant(0), left)
    bc_R = DirichletBC(V, Constant(0), right)
    bc_T = DirichletBC(V, Constant(0), top)
    bc_B = DirichletBC(V, Constant(0), bottom)

    ## Apply boundary conditions
    bcs = [bc_L,bc_R,bc_T,bc_B]


    for bc in bcs:
        bc.apply(A_g)

    A = A_g


    return A

################################################################################

####              Stochastic Part of Code                        #

################################################################################


### Ajit's Part Used
## Spatialy varying KLE/PCE modes
class MyExpression(Expression):
    def __init__(self, params, **kwargs):
        self.index  = params[0]
        self.ndim   = params[1]
        self.sIndex = params[2].astype(int)
        self.mIndex = params[3].astype(int)
        self.PCE_A  = params[4]
        self.xi     = params[5]
        self.PCE_u  = params[6]
        self.PC1D   = params[7]
        self.meanl  = params[8]
        self.sigma  = params[9]

    def eval(self, values, x):
        ## For bx=by=1 good for upto 15 RVs cases (unit square)
        multipliers = [0.92184, 0.49248, 0.29374, 0.20437, 0.15576]
        omegas = [1.30654, 3.67319, 6.58462, 9.63168, 12.72324]

        # Log normal Mean and Standard Deviation of Underlaying Gaussian
        # meanl = self.meanl     # similar to meanc in Fortran package
        # sigma = self.sigma      # similar to sigma in Fortran package

### Trunctaed PCE of lognormal process
## Automation to n-RVs : ----------------------------------------------------------------
## KLE: obtain KLE terms
        g = []
        for i in range(self.ndim):
            Xindex = self.sIndex[i,0]
            # print('Xindex', Xindex)
            ### First Dimesnion
            if (Xindex % 2) == 0:
                ##print("even")
                gg1 = multipliers[Xindex-1] * (sin(omegas[Xindex-1]*(x[0]-0.5)))
            else:
                ##print("odd")
                gg1 = multipliers[Xindex-1] * (cos(omegas[Xindex-1]*(x[0]-0.5)))

            ### Second Dimension
            Yindex = self.sIndex[i,1]
            if (Yindex % 2) == 0:
                ##print("even")
                gg2 = multipliers[Yindex-1] * (sin(omegas[Yindex-1]*(x[1]-0.5)))
            else:
                ##print("odd")
                gg2 = multipliers[Yindex-1] * (cos(omegas[Yindex-1]*(x[1]-0.5)))

            g.append(self.sigma*gg1*gg2)


## PCE: obtain PCE terms
        kappa_coeff = []
        for ipce in range(PCE_A):
            newY = 1.0
            # i = self.index
            for j in range(self.ndim):
                idx = mIndex[ipce,j]
                if idx == 0:
                    yy = 1.0
                elif idx == 1:
                    yy = g[j]
                else:
                    varFacto = nfactorial(idx)
                    yy = (g[j]**(idx))/np.sqrt(varFacto)

                ## Keep all indentations properly (there was error here due to indentation)
                newY =  newY*yy

            # print('kappa_coeff is', kappa_coeff)
            # print('ipce2 is', ipce)
            kappa_coeff.append(self.meanl * newY)



        lterms = np.zeros([self.PCE_u])

        lterms = multiPC(self.ndim,self.PCE_u, mIndex, PC1D, xi)
        #### Multiplying chaos terms with the Gaussian Sample for MCS


        # xi_1 = self.xi[0]
        # xi_2 = self.xi[1]
        # xi_3 = self.xi[2]
        # # print('xi1 and xi2 is', xi_1,xi_2)

        # lterms = np.zeros([self.PCE_u])

        # if self.ndim == 2:


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


        #####
        # print(lterms)

        # exit()
        kappa = 0

        for n_l in range(PCE_A):

            kappa = kappa + kappa_coeff[n_l]*lterms[n_l]


        values[0] = kappa

##---------------------------------------------------------------------------
## Stochastic variational form
def stoVartional(params, u1, w, f):
    kleID = params[0]
    # print('KLE id', kleID)
    cd = MyExpression(params, degree=0)
    # print('cd is',cd)

    m = u1*w*dx
    a = cd* inner(nabla_grad(u1),nabla_grad(w))*dx
    l = f*w*dx
    return a,m,l


def newmark_update(u, u0, v0, a0, beta, gamma, dt):

    u_vec, u0_vec = u.vector(), u0.vector()
    v0_vec, a0_vec = v0.vector(), a0.vector()

 # Update acceleration and velocity
    a_vec = (1.0/(2.0*beta))*( (u_vec - u0_vec -
    v0_vec*dt)/(0.5*dt*dt) - (1.0-2.0*beta)*a0_vec )

    v_vec = dt*((1.0-gamma)*a0_vec + gamma*a_vec) + v0_vec
    v0.vector()[:], a0.vector()[:] = v_vec, a_vec
    u0.vector()[:] = u.vector()



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
#########################################################################


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

###############################    Main Code Starts ##################################################################




if(ip ==0):
    t1 = time.time()


T = 3;
dt = 6.5e-3;
beta = 0.25;
gamma = 0.5;

t = 0.0
nt = np.int(T/dt)


# Reading Mesh information
# mpath = "./meshData/square_acoustic.xml"
# mesh = Mesh(mpath)

mesh = Mesh(mpi_comm_self())
hdf = HDF5File(mesh.mpi_comm(), "./../meshData/file.h5", "r")
hdf.read(mesh, "/mesh", False)
cd = CellFunction("size_t", mesh)
hdf.read(cd, "/cd")
fd = FacetFunction("size_t", mesh)
hdf.read(fd, "/fd")

ndof = mesh.num_vertices()

if(ip ==0):
    print('number of degrees of freedom',ndof)



# ndim = 2
# mean_g = 0
# sigma_g = 0.1


pcedatapath = "./../klePceData/pcedata.dat"
pcedata = np.genfromtxt(pcedatapath)
order_u = pcedata[0].astype(int)     ## order of output PCE (input==2)
dim_L   = pcedata[1].astype(int)     ## number of KLE modes (RVs)
PCE_u   = pcedata[2].astype(int)     ## number of PCE outputs
PCE_A   = pcedata[3].astype(int)    ## number of PCE inputs

### Log normal mean exp(mean_g + 0.5*sigma_g^2)
meanl = 1.046

sigma_g = 0.3


if(ip ==0):
    print('Random Dimension', dim_L)
    print('number of PC input', PCE_A)
    print('Log normal mean',meanl)
    print('Gaussian sd',sigma_g)


## Load KLE/PCE Indices for Automation
sIndexPath = "./../klePceData/sortIndex.dat"
sIndex = np.genfromtxt(sIndexPath).astype(int)

mIndexPath = "./../klePceData/multiIndex.dat"
mIndex = np.genfromtxt(mIndexPath).astype(int)

# print('sort index is', sIndex)


def nfactorial(nf):
    if nf == 0:
        return 1
    else:
        return nf * nfactorial(nf-1)


if(ip ==0):
    nSamples = 10000
    print('number of samples',nSamples)
    NP = 400
    chunks = int(nSamples/NP)
    print('chunks is ', chunks)
else:
    chunks = None
    NP = None

chunks = comm.bcast(chunks, root = 0)
NP = comm.bcast(NP, root = 0)


u_MCS = np.zeros([ndof,nt+1,chunks])
u_sum = np.zeros([ndof,nt+1])
u_MCS_sd = np.zeros([ndof,nt+1])

u_pdf = np.zeros([9,chunks])

t_s1 = time.time()


for k in range(chunks):


        # Function space definition over the mesh
    V = FunctionSpace(mesh, "CG", 1)   #Continuous Galerkin for the displacement field

    # Test and Trial function
    u1, w = TrialFunction(V), TestFunction(V)

    # Initialization of fields (displacement, velocity, acceleration)
    u0, v0, a0 = Function(V), Function(V), Function(V)


      #Initial conditions
    f  = Constant((0.0))

    ui  = Expression(("1*exp(-100*(pow(x[0]-0.7,2)+pow(x[1]-0.7,2)))"),degree=1)


    ##### Initial condition ######

    # ui  = Expression(("sin(m*pi*x[0])*sin(n*pi*x[1])"),degree=1,m = 2,n = 1)

    u0 = interpolate(ui, V)

    v0 = interpolate(Constant(0.0), V)

    a0 = interpolate(Constant(0.0), V)

    c = 1 # wave velocity


######################################

    xi = npr.randn(dim_L)

    print('Process id is',ip)
    print('xi is',xi)

    PC1D = OnedPC(dim_L, xi)

    params = ([k, dim_L, sIndex, mIndex, PCE_A, xi, PCE_u, PC1D,meanl,sigma_g])

    ## Invoke variation fomrulaiton for each sample
    a,m,l = stoVartional(params, u1, w, f)


##################
    # m = u1*w*dx
    # a = c* c * inner(nabla_grad(u1),nabla_grad(w))*dx
    # l = f*w*dx
##################


    if(k ==0):
    ## Invoke Deterministic Assembly procedure for each sample
        A,M,b = detAssembly(a,m,l)     ## New Assembly
    else:
        A = detAssemblyA(a)


    # t1 = time.time()
    # tt = t1-t0
    # ### Damped
    #C = 0.0 * M + 0.00 * A;
    # print('time taken for 1 sample assembling',tt )

    C = 0.5445 * M + 0.0174 * A;
    #


    K_t = (M/(beta*dt*dt)) + (C *gamma/(beta*dt)) +  A

    # print(K_t)

    u = Function(V)

    tseries = np.zeros([ndof,nt+1])

    count = 0

    while count < nt+1:

       t += dt

       # print('t is', t)
       # print('b is', b)


       f_m =  ((u0.vector() + dt*v0.vector()) /(beta*dt*dt)) + ((1.0-2.0*beta)*a0.vector()/(2*beta))
       f_c = gamma * dt * f_m - v0.vector() - (1-gamma)*a0.vector()*dt
       F_t = b + M * f_m + C * f_c
       #F_t = b + M * ( ((u0.vector() + dt*v0.vector()) /(beta*dt*dt)) + ((1.0-2.0*beta)*a0.vector()/(2*beta)) )

       solve(K_t,u.vector(),F_t)

       newmark_update(u, u0, v0, a0, beta, gamma, dt)

       # file << u
       tseries[0:ndof,count]= u.compute_vertex_values()

       count+=1
       # print(count)


    u_pdf[0,k] = tseries[400,149]
    u_pdf[1,k] = tseries[400,299]
    u_pdf[2,k] = tseries[400,461]

    u_pdf[3,k] = tseries[401,149]
    u_pdf[4,k] = tseries[401,299]
    u_pdf[5,k] = tseries[401,461]

    u_pdf[6,k] = tseries[402,149]
    u_pdf[7,k] = tseries[402,299]
    u_pdf[8,k] = tseries[402,461]


#### Inside sample loop
    u_MCS[:,:,k] = tseries


    ###### Summing up Monte Carlo Solutions to find Mean solution
    u_sum = u_sum + tseries

################

spath = './u_pdf_'+ str(ip) +'.mat'
sio.savemat(spath, {'u_pdf':u_pdf})

# spath = './u_MCS.mat'
# sio.savemat(spath, {'u_MCS':u_MCS})

# spath = './u_mean.mat'
# sio.savemat(spath, {'u_mean':u_mean})

u_MCS_mean = u_sum/chunks


for j in range(chunks):

    u_MCS_sd = u_MCS_sd + np.power((u_MCS[:,:,j] - u_MCS_mean),2)


u_MCS_sd = u_MCS_sd/(chunks-1)


mu_sum = np.zeros([ndof,nt+1])
sd_sum = np.zeros([ndof,nt+1])

mu_sum = comm.allreduce(u_MCS_mean, op=MPI.SUM)

comm.Reduce(u_MCS_sd, sd_sum, op=MPI.SUM, root=0)



mu_exact = mu_sum/NP

var_mu = np.power((u_MCS_mean - mu_exact),2)

var_mu_sum = np.zeros([ndof,nt+1])

comm.Reduce(var_mu, var_mu_sum, op=MPI.SUM, root=0)



if(ip == 0):
    sd_sum = sd_sum * (chunks-1)/(nSamples-1)

    v_mu = var_mu_sum * (chunks)/(nSamples-1)

    v_exact = sd_sum + v_mu

    sd_exact = np.sqrt(v_exact)
    # print('mu sum is', mu_sum)
    t2 = time.time()
    print('time taken is', t2-t1)
    print("....................................")
    print("................Success.............")
    print(".....................................")

    spath = './nS.mat'
    sio.savemat(spath, {'nS':nSamples})


    spath = './u_MCS_mean.mat'
    sio.savemat(spath, {'u_MCS_mean':mu_exact})


    spath = './u_MCS_sd.mat'
    sio.savemat(spath, {'u_MCS_sd':sd_exact})















# u_MCS_sd = np.sqrt(u_MCS_sd)

# spath = './nS.mat'
# sio.savemat(spath, {'nS':nSamples})


# spath = './u_MCS_mean.mat'
# sio.savemat(spath, {'u_MCS_mean':u_MCS_mean})


# spath = './u_MCS_sd.mat'
# sio.savemat(spath, {'u_MCS_sd':u_MCS_sd})


















