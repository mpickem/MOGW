F2PY = f2py
FC = ifort
FFLAGS += -O3 -g -lmkl_rt

all: comp.so glb.so ham.so index_reference.o lapack_module.o

comp.so: computation_module.F90 global.o index_reference.o hamiltonian_module.o lapack_module.o
	$(F2PY) -c -m comp $^ --fcompiler=intelem --quiet -lmkl_rt

glb.so: global.f90
	$(F2PY) -c -m glb $^ --fcompiler=intelem --quiet

ham.so: hamiltonian_module.f90 global.f90
	$(F2PY) -c -m ham $^ --fcompiler=intelem --quiet

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS) -fPIC

index_reference.o: global.o
hamiltonian_module: global.o
lapack_module.o: global.o


.PHONY: clean
clean:
	rm -rf *.so *.mod *.o
