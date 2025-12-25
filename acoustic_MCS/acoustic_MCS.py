


### Acoustic Wave Propagation Problem - Serial ######

##### Monte Carlo Simulation using Fortran Decomposed Mesh and Fenics Assembly

#### Copyright (C) Sudhi P V #####

### June -30 - 2020





from dolfin import *
import numpy as np
import time
import os as os
import scipy.io as sio
import numpy.random as npr


## Initialize MPI for python
# comm = MPI.COMM_WORLD
# ip = comm.Get_rank()
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
        self.PCE_u     = params[6]

    def eval(self, values, x):
        ## For bx=by=1 good for upto 15 RVs cases (unit square)
        multipliers = [0.92184, 0.49248, 0.29374, 0.20437, 0.15576]
        omegas = [1.30654, 3.67319, 6.58462, 9.63168, 12.72324]

        # Log normal Mean and Standard Deviation of Underlaying Gaussian
        meanl = 1.005     # similar to meanc in Fortran package
        sigma = 0.1      # similar to sigma in Fortran package

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

            g.append(sigma*gg1*gg2)


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
            kappa_coeff.append(meanl * newY)


        #### Multiplying chaos terms with the Gaussian Sample for MCS
        xi_1 = self.xi[0]
        xi_2 = self.xi[1]
        xi_3 = self.xi[2]
        # print('xi1 and xi2 is', xi_1,xi_2)

        lterms = np.zeros([self.PCE_u])

        if self.ndim == 2:


            lterms[0] = 1.0
            lterms[1] = xi_1
            lterms[2] = xi_2
            lterms[3] = pow(xi_1,2) - 1;
            lterms[4] = xi_1*xi_2;
            lterms[5] = pow(xi_2,2) - 1;
            lterms[6] = pow(xi_1,3) - 3*xi_1;
            lterms[7] = xi_1*xi_1*xi_2 - xi_2;
            lterms[8] = xi_2*xi_2*xi_1 - xi_1;
            lterms[9] = pow(xi_2,3) - 3*xi_2;


        if self.ndim == 3:

            lterms[0] = 1.0
            lterms[1] = xi_1;
            lterms[2] = xi_2;
            lterms[3] = xi_3;
            lterms[4] = pow(xi_1,2) - 1;
            lterms[5] = xi_1*xi_2;
            lterms[6] = xi_1*xi_3;
            lterms[7] = pow(xi_2,2) - 1;
            lterms[8] = xi_2*xi_3;
            lterms[9] = pow(xi_3,2) - 1;
            lterms[10] = pow(xi_1,3) - 3*xi_1;
            lterms[11] = pow(xi_1,2)*xi_2 - xi_2;
            lterms[12] = pow(xi_1,2)*xi_3 - xi_3;
            lterms[13] = pow(xi_2,2)*xi_1 - xi_1;
            lterms[14] = xi_1*xi_2*xi_3;
            lterms[15] = pow(xi_3,2)*xi_1 - xi_1;
            lterms[16] = pow(xi_2,3) - 3*xi_2;
            lterms[17] = pow(xi_2,2)*xi_3 - xi_3;
            lterms[18] = pow(xi_3,2)*xi_2 - xi_2;
            lterms[19] = pow(xi_3,3) - 3*xi_3;


        #####

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



###############################    Main Code Starts ##################################################################


T = 0.5;
dt = 2.5e-3;
beta = 0.25;
gamma = 0.5;

t = 0.0
nt = np.int(T/dt)


# Reading Mesh information
# mpath = "./meshData/square_acoustic.xml"
# mesh = Mesh(mpath)

mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "./meshData/file.h5", "r")
hdf.read(mesh, "/mesh", False)
cd = CellFunction("size_t", mesh)
hdf.read(cd, "/cd")
fd = FacetFunction("size_t", mesh)
hdf.read(fd, "/fd")

ndof = mesh.num_vertices()

print('number of degrees of freedom',ndof)



# ndim = 2
# mean_g = 0
# sigma_g = 0.1


pcedatapath = "./klePceData/pcedata.dat"
pcedata = np.genfromtxt(pcedatapath)
order_u = pcedata[0].astype(int)     ## order of output PCE (input==2)
dim_L   = pcedata[1].astype(int)     ## number of KLE modes (RVs)
PCE_u   = pcedata[2].astype(int)     ## number of PCE outputs
PCE_A   = pcedata[3].astype(int)    ## number of PCE inputs

## Load KLE/PCE Indices for Automation
sIndexPath = "./klePceData/sortIndex.dat"
sIndex = np.genfromtxt(sIndexPath)

mIndexPath = "./klePceData/multiIndex.dat"
mIndex = np.genfromtxt(mIndexPath)

# print('sort index is', sIndex)


def nfactorial(nf):
    if nf == 0:
        return 1
    else:
        return nf * nfactorial(nf-1)


nSamples = 2

print('number of samples',nSamples)

u_MCS = np.zeros([ndof,nt+1,nSamples])
u_mean = np.zeros([ndof,nt+1])
u_MCS_sd = np.zeros([ndof,nt+1])

t_s1 = time.time()



for k in range(nSamples):


        # Function space definition over the mesh
    V = FunctionSpace(mesh, "CG", 1)   #Continuous Galerkin for the displacement field

    # Test and Trial function
    u1, w = TrialFunction(V), TestFunction(V)

    # Initialization of fields (displacement, velocity, acceleration)
    u0, v0, a0 = Function(V), Function(V), Function(V)


    def dot_u(u):
      return (gamma/(beta*dt))*(u - u0) - (gamma/beta - 1.0)*v0 - dt*(gamma/(2.0*beta) - 1.0)*a0

    def ddot_u(u):
      return (1.0/(beta*dt**2))*(u - u0 - dt*v0) - (1.0/(2.0*beta) - 1.0)*a0


      #Initial conditions
    f  = Constant((0.0))

    ui  = Expression(("1*exp(-100*(pow(x[0]-0.5,2)+pow(x[1]-0.5,2)))"),degree=1)


    ##### Initial condition ######

    # ui  = Expression(("sin(m*pi*x[0])*sin(n*pi*x[1])"),degree=1,m = 2,n = 1)

    u0 = interpolate(ui, V)

    v0 = interpolate(Constant(0.0), V)

    a0 = interpolate(Constant(0.0), V)

    c = 1 # wave velocity



######################################

    xi = npr.randn(dim_L)

    print(xi)

    params = ([k, dim_L, sIndex, mIndex, PCE_A, xi, PCE_u])

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
    C = 0.0 * M + 0.00 * A;
    # print('time taken for 1 sample assembling',tt )

    # C = 0.8624 * M + 0.0013 * A;
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

       F_t = b + M * ( ((u0.vector() + dt*v0.vector()) /(beta*dt*dt)) + ((1.0-2.0*beta)*a0.vector()/(2*beta)) )

       solve(K_t,u.vector(),F_t)

       newmark_update(u, u0, v0, a0, beta, gamma, dt)

       # file << u
       tseries[0:ndof,count]= u.compute_vertex_values()

       count+=1
       # print(count)



#### Inside sample loop
    u_MCS[:,:,k] = tseries


    ###### Summing up Monte Carlo Solutions to find Mean solution
    u_mean = u_mean + tseries

################


spath = './u_MCS.mat'
sio.savemat(spath, {'u_MCS':u_MCS})

# spath = './u_mean.mat'
# sio.savemat(spath, {'u_mean':u_mean})


t_s2 = time.time()

total = t_s2 - t_s1

print('total time taken',total)

u_MCS_mean = u_mean/nSamples


for j in range(nSamples):

    u_MCS_sd = u_MCS_sd + np.power((u_MCS[:,:,j] - u_MCS_mean),2)


u_MCS_sd = u_MCS_sd/(nSamples-1)

u_MCS_sd = np.sqrt(u_MCS_sd)


spath = './nS.mat'
sio.savemat(spath, {'nS':nSamples})


spath = './u_MCS_mean.mat'
sio.savemat(spath, {'u_MCS_mean':u_MCS_mean})


spath = './u_MCS_sd.mat'
sio.savemat(spath, {'u_MCS_sd':u_MCS_sd})


















