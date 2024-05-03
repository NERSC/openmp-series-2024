#!/bin/bash
#SBATCH --time=05:00
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --constraint=cpu
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=256
#SBATCH --job-name=omp_for_scaling
#SBATCH --output=omp_for_scaling.%j.txt

# ================ optional OMP - settings (thread pinning)
# thread-placement options: cores, threads, sockets
# export OMP_PLACES=cores 
# echo Thread Placement: ${OMP_PLACES}

# thread-binding options: close, spread
# export OMP_PROC_BIND=close
# echo Thread Binding Policy: ${OMP_PROC_BIND}

# =============== Run with different no of threads
NTHREADS="1 2 3 4 6 8 12 16 24 32 40 48"
#NTHREADS="1 2 3 4 6 8 12"
for i in $NTHREADS
do
    echo Number of Threads = $i
    OMP_NUM_THREADS=$i make -s run
    echo
done