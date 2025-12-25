#!/bin/bash
#SBATCH --account=def-asarkar
#SBATCH --time=0-2:30
#SBATCH --mem-per-cpu=4500M
#SBATCH --nodes=10
#SBATCH --tasks-per-node=40
#SBATCH --job-name=slurm
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=sudhipv@cmail.carleton.ca
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


source fenics_activate.sh

num_sim="10"

for i in $(seq 1 $num_sim);do

mkdir ./3rv_sigma3_pdf_$i

cp ./acoustic_pdf.py ./3rv_sigma3_pdf_$i

cd ./3rv_sigma3_pdf_$i

mpiexec -n 400 python3 acoustic_pdf.py
#mpiexec -n 200 python3 acoustic_parallel_saving_MCS.py

done


# mpiexec -n 400 python3 acoustic_parallel_MCS.py
#mpiexec -n 200 python3 acoustic_parallel_saving_MCS.py








