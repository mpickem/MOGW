#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
from mpi4py import MPI
import h5py
import sys
import os
import time
import matplotlib.pyplot as plt

from src.glb import mglobal
from src.ham import mhamil
from src.comp import mcomp

'''
Main wrapper for the underlying fortran files
found in src. Calculation of the Selfenergy approximation GW
in the Wannier Basis
'''

# MPI.Init()
# mpi communicator, rank, size + data distribution
comm = MPI.COMM_WORLD
master = 0 # for easy access
mpi_rank = comm.Get_rank()
mpi_size = comm.Get_size()


# reading the Hamiltonian stuff
pwd = os.getcwd()
hmlt = pwd + '/input_srvo3_full/HMLT'
hmlt_kpq = pwd + '/input_srvo3_full/HMLT.index.kpq'
hmlt_mq = pwd + '/input_srvo3_full/HMLT.index.mq'

for i,j in enumerate(hmlt):
  mhamil.file_hmlt[i] = j
for i,j in enumerate(hmlt_kpq):
  mhamil.file_hmlt_kpq[i] = j
for i,j in enumerate(hmlt_mq):
  mhamil.file_hmlt_mq[i] = j

mhamil.read_hamiltonian()
mhamil.read_bzindices()

data_interval = mhamil.nkp # we get this number from the hamiltonian and/or input
ndim = mhamil.ndim
nkp = mhamil.nkp
mcomp.wtkp = mhamil.wtkp

nw = 50
mcomp.beta = 10

# MPI distribution
displ=[] # displacement
displ.append(0)
rct=[]   # receive count

# distribution of data_interval
for i in xrange(mpi_size-1):
  rct.append((data_interval-displ[i])//(mpi_size-i))
  displ.append(rct[i]+displ[i])
rct.append(data_interval - displ[mpi_size-1])

if (mpi_rank == master):
  print('rct: ', rct)
  print('displ: ', displ)
  print()
  sys.stdout.flush()

rct = np.array(rct)
displ = np.array(displ)

ikstart = displ[mpi_rank] + 1
ikend   = displ[mpi_rank] + rct[mpi_rank]

comm.barrier()

print(ikend,ikstart)

Giw = np.zeros((ndim,ndim,2*nw,nkp), dtype=np.complex128, order='F')
Selfenergy = np.zeros((ndim,ndim,2*nw,nkp), dtype=np.complex128, order='F')

mu = 0.0
mcomp.compute_giw(Giw,mhamil.h,Selfenergy,mu,ikstart,ikend)


Giw_gather = np.zeros_like(Giw, dtype=np.complex128, order='F')
comm.Allgatherv([Giw, rct[mpi_rank]*ndim**2*2*nw, displ[mpi_rank]*ndim**2*2*nw, MPI.COMPLEX16], \
                [Giw_gather,rct,displ, MPI.COMPLEX16])

if (mpi_rank == master):
  plt.plot(Giw_gather[0,0,:,0].real)
  plt.plot(Giw_gather[0,0,:,0].imag)
  plt.show()


nmatrix = np.zeros((ndim,ndim), dtype=np.complex128, order='F')
occ = mcomp.compute_n(Giw_gather[:,:,:,:],nmatrix)
print(occ.real)


# a = mhamil.h.reshape((3,3,20,20,20)).real # this is fine since we have an F ordered array
