#!/bin/bash
#SBATCH -p defq
#SBATCH -n 60
#SBATCH --ntasks-per-node=60
#SBATCH -J ctint
#SBATCH -o stdout.%J
#SBATCH -e stderr.%J

module load openmpi/3.1.5/gcc-9.3.0

python param.py
python gen_G0.py
date
# Use "mpirun" command in your system. In this case, we use 60 processes.
mpirun -np 60 $HOME/opt/CT-INT/bin/ctint_real input.ini > output
date

# read data from input.out.h5
python read_G.py

# plot data
python plot.py
