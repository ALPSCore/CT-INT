#!/bin/sh
#BSUB -n 120
#BSUB -o debug
#BSUB -J debug

python param.py
python gen_G0.py
date
# Use "mpirun" command in your system. In this case, we use 120 processes.
/usr/share/lava/1.0/linux2.6-glibc2.12-x86_64/bin/intelmpi-mpirun -np 120 /home/shinaoka/build/ct-int-alpscore/ctint_real input.ini > output
date

# read data from input.out.h5
python read_G.py

# plot data
python plot.py
