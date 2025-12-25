from dolfin import *

mesh = Mesh("square_acoustic.xml");
cd=MeshFunction('size_t',mesh,"square_acoustic_physical_region.xml");
fd=MeshFunction('size_t',mesh,"square_acoustic_facet_region.xml");
hdf = HDF5File(mesh.mpi_comm(), "file.h5", "w")
hdf.write(mesh, "/mesh")
hdf.write(cd, "/cd")
hdf.write(fd, "/fd")
