## KLE/PCE data generator and importer
## This take one command line input after number of dimensions (RVs)
cd klePceData/
gfortran KLE_PCE_Data_commandLineArg.F90
./a.out
cd -

## GMSH data generator
## Change to respective folder

cd meshData/
cp foo.geo square_acoustic.geo


## Adjust the mesh density parameter "lc"
## lc=old;/lc=new for foo*.geo

sed -i -e 's/lc=0.02;/lc=0.01;/g' square_acoustic.geo

## select the number of partitions

## create mesh file::
gmsh -2 square_acoustic.geo -o square_acoustic.msh


cd -


