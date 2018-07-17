#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
from mpi4py import MPI
import h5py
import sys
import os
import time

from src.glb import mglobal
from src.ham import mhamil
from src.comp import mcomp

'''
Main wrapper for the underlying fortran files
found in src. Calculation of the Selfenergy approximation GW
in the Wannier Basis
'''

# mpi communicator, rank, size + data distribution
comm = MPI.COMM_WORLD
master = 0 # for easy access
mpi_rank = comm.Get_rank()
mpi_size = comm.Get_size()

data_interval = 31 # we get this number from the hamiltonian and/or input

displ=[] # displacement
displ.append(0)
rct=[]   # receive count

# distribution of data_interval
for i in xrange(mpi_size-1):
  rct.append((data_interval-displ[i])//(mpi_size-i))
  displ.append(rct[i]+displ[i])
rct.append(data_interval - displ[mpi_size-1])

rct = np.array(rct)
displ = np.array(displ)

mglobal.ikstart = displ[mpi_rank] + 1   # fortran
mglobal.ikend   = displ[mpi_rank] + rct[mpi_rank]

pwd = os.getcwd()
hmlt = pwd + '/input_srvo3_full/HMLT'
hmlt_kpq = pwd + '/input_srvo3_full/HMLT.index.kpq'
hmlt_mq = pwd + '/input_srvo3_full/HMLT.index.mq'

mglobal.clear_file_arrays()

for i,j in enumerate(hmlt):
  mglobal.file_hmlt[i] = j
for i,j in enumerate(hmlt_kpq):
  mglobal.file_hmlt_kpq[i] = j
for i,j in enumerate(hmlt_mq):
  mglobal.file_hmlt_mq[i] = j

print(mglobal.file_hmlt)

mglobal.giw = np.ones((2,2,2,2), dtype=np.complex128)

mhamil.read_hamiltonian()

sys.exit()

mhamil.read_bzindices()

if (mpi_rank == master):
  print(mglobal.h)
