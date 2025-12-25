#!/bin/bash
#SBATCH --account=def-asarkar
#SBATCH --time=0-00:40
#SBATCH --mem-per-cpu=3700M
#SBATCH --nodes=4
#SBATCH --tasks-per-node=48
#SBATCH --job-name=slurm
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=sudhipv@cmail.carleton.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


source fenics_activate.sh

#mpiexec -n 4 python3 nisp_parallel1.py

#mpiexec -n 3 python3 nisp_parallel2.py

#rm u_nisp_part1_mean.mat

#python3 nisp_serial.py
mpiexec -n 159 python3 nisp_SparseQuad.py
#mpiexec -n 4 python3 nisp_Sparse_NEW.py
#mpiexec -n 200 python3 acoustic_parallel_saving_MCS.py









