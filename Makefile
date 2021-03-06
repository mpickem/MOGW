#Makefile for VSC3

#Fortran Compiler
FC        = mpiifort

# FFLAGS = -g -O0 -heap-arrays 10 -traceback -check bounds -check uninit -fpe-all=3 # ifort debug
#Fortran Flags
FFLAGS    = -O3 -g -fpp -DMPI  # ifort production

#Library Directories
LDIR      = -L/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.12/intel-14.0.2/lib/ \
						-L/cm/shared/apps/intel/composer_xe_2015.2.164/mkl/lib/intel64

#Libraries
LIBS     += -lmkl_rt -lhdf5_fortran -lhdf5hl_fortran

#Include Directories
IDIR     += -I/opt/sw/x86_64/glibc-2.12/ivybridge-ep/hdf5/1.8.12/intel-14.0.2/include

#Fortran Sources
FSOURCES  = aux.f90 computation_functions.f90 four.f90 gw.f90 \
					 hamiltonian_module.f90 hdf5_module.f90 index_reference.f90 \
					 io.f90 lapack_module.f90 mpi_org.f90 read_functions.f90 \
					 vq_module.f90

#Pattern Substitution from *.f90 to *.o
OBJ = $(FSOURCES:.f90=.o)

#Compilation Rules
.PHONY: all
all: gw

gw: $(OBJ)
	$(FC) $^ -o $@ $(FFLAGS) $(LDIR) $(LIBS) $(IDIR)

gw.o:                    aux.o computation_functions.o four.o hamiltonian_module.o hdf5_module.o index_reference.o io.o \
												 lapack_module.o mpi_org.o read_functions.o vq_module.o
hdf5_module.o:           aux.o hamiltonian_module.o
vq_module.o:             aux.o hamiltonian_module.o hdf5_module.o
index_reference.o:       aux.o hamiltonian_module.o
computation_functions.o: aux.o index_reference.o hamiltonian_module.o lapack_module.o mpi_org.o
read_functions.o:        aux.o index_reference.o hamiltonian_module.o mpi_org.o hdf5_module.o vq_module.o
io.o:                    aux.o hamiltonian_module.o
hamiltonian_module.o:    aux.o
lapack_module.o:         aux.o
mpi_org.o:               aux.o

$(OBJ): %.o : %.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(IDIR)

.PHONY: clean
clean :
	rm -f gw *.mod *.o
