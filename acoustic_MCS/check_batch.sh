

#source fenics_activate.sh

num_sim="2"

for i in $(seq 1 $num_sim);do

mkdir ./3rv_sigma3_pdf_$i

#mpiexec -n 400 python3 acoustic_parallel_MCS.py
#mpiexec -n 200 python3 acoustic_parallel_saving_MCS.py

done










