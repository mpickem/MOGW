#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
from mpi4py import MPI
import h5py
import sys
import os
import time

from src.comp import mcomp

'''
Main wrapper for the underlying fortran files
found in src. Calculation of the Selfenergy approximation GW
in the Wannier Basis
'''

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
  mcomp.file_hmlt[i] = j
for i,j in enumerate(hmlt_kpq):
  mcomp.file_hmlt_kpq[i] = j
for i,j in enumerate(hmlt_mq):
  mcomp.file_hmlt_mq[i] = j

mcomp.read_hamiltonian()
mcomp.read_bzindices()

data_interval = mcomp.nkp
ndim = mcomp.ndim
nkp = mcomp.nkp

mcomp.nw = nw = 50
mcomp.beta = beta = 10

# MPI distribution
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

if (mpi_rank == master):
  print('rct: ', rct)
  print('displ: ', displ)
  print()
  sys.stdout.flush()

mcomp.ikstart = displ[mpi_rank] + 1
mcomp.ikend = displ[mpi_rank] + rct[mpi_rank]

# Giw = np.zeros((ndim,ndim,2*nw,nkp), dtype=np.complex128, order='F')
# Selfenergy = np.zeros((ndim,ndim,2*nw,nkp), dtype=np.complex128, order='F')

mcomp.gkiw = np.zeros((ndim,ndim,2*nw,nkp), dtype=np.complex128, order='F')
mcomp.skiw = np.zeros((ndim,ndim,2*nw,nkp), dtype=np.complex128, order='F')

mcomp.compute_gkiw(0.0)

gkiw_gather = np.zeros_like(mcomp.gkiw, dtype=np.complex128, order='F')

mpi_ctemp = np.zeros((ndim,ndim,nkp), dtype=np.complex128, order='F')
for i in xrange(2*nw):
  a = np.asfortranarray(mcomp.gkiw[:,:,i,:][()])
  comm.Allgatherv([a, rct[mpi_rank]*ndim**2, displ[mpi_rank]*ndim**2, MPI.COMPLEX16], \
                  [mpi_ctemp,rct*ndim**2,displ*ndim**2,MPI.COMPLEX16])
  gkiw_gather[:,:,i,:] = mpi_ctemp[()]
del mpi_ctemp

mcomp.gkiw = gkiw_gather

occ = mcomp.compute_n()
print(occ.real)


# a = mhamil.h.reshape((3,3,20,20,20)).real # this is fine since we have an F ordered array
