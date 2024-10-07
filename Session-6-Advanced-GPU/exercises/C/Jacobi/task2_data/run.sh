#!/bin/bash

### Request GPU node
#SBATCH -N 1
#SBATCH -C gpu

### specify NERSC project: ntrain5 or your own project such as m1234
#SBATCH -A ntrain5

### Request 1 GPU and 32 CPU cores on the GPU node (share with others)
#SBATCH --qos=shared
#SBATCH -G 1
#SBATCH -c 32
 
### Job name
#SBATCH -J Task1

### File / path where STDOUT & STDERR will be written
#SBATCH -o Task1.%J.out
#SBATCH -e Task1.%J.err
 
### Request the time you need for execution in minutes
### The format for the parameter is: hour:minute:seconds
#SBATCH -t 00:10:00

### setup environment
module load PrgEnv-nvidia

### Compile & Execute your application
make clean
make
#make run
OMP_NUM_THREADS=1 ./jacobi.gpu
#OMP_NUM_THREADS=1 ./jacobi.sol.gpu
