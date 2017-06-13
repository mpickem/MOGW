# Makefile for VSC3

FORTRAN = mpiifort 

# FFLAGS = -g -O0 -heap-arrays 10 -traceback -check bounds -check uninit -fpe-all=3 # ifort debug
FFLAGS = -O3 -g -fpp -DMPI  # ifort production

LIBS = -L$(LIBRARY_PATH) -lmkl_rt -lhdf5_fortran -lhdf5hl_fortran

# -L$(MKLROOT)/lib/intel64
#-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread 
#-lmkl_lapack -lguide

PROG = gw_dmft
OBJS = aux.o lapack_module.o hamiltonian_module.o four.o mpi_org.o hdf5_module.o vq_module.o gw.o io.o
# parameters_module.o 

all : $(PROG) $(OBJS)

$(PROG) : $(OBJS)
	$(FORTRAN) -o $(PROG) $(FFLAGS) $(OBJ) $(LDFLAGS) $(LIBS) $^

%.o : %.F90
	$(FORTRAN) -c $(FFLAGS) $< 

clean :
	rm -f $(PROG) *.mod *.o
