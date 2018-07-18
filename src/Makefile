F2PY = f2py
FC = ifort
FFLAGS += -O3 -g -lmkl_rt

all: comp.so glb.so lapack_module.o

comp.so: computation_module.F90 global.o lapack_module.o
	$(F2PY) -c -m comp $^ --fcompiler=intelem --quiet -lmkl_rt

glb.so: global.f90
	$(F2PY) -c -m glb $^ --fcompiler=intelem --quiet

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS) -fPIC

lapack_module.o: global.o

.PHONY: clean
clean:
	rm -rf *.so *.mod *.o
