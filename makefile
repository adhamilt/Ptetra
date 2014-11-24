OBJS=   ptetra.o
FLAGS=  -O -r8 -warn
FLAGS=  -O0 -g -C -r8 -check all -traceback
FLAGS=  -O3 -r8 -warn
FLAGS=  -O3 -fdefault-real-8 -fdefault-double-8
#FLAGS=  -O0 -fdefault-real-8 -fdefault-double-8 -fcheck=all

SPARSKIT= -L/home/richard/lib/ -lskit
BLAS= -L/home/richard/lib/ -lblas
#BLAS=  -L/opt/intel/mkl/10.0.1.014/lib/32/ -lmkl_core.so -lmkl_blas95.a

#put the libblas.a and libskit files directly in this directory or somehow
#arrange for the linker to know where to find them

#F90=   g95
F90=    ifort
F90=    gfortran

ptetra: ptetra.o
	$(F90) -o ptetra $(FLAGS) $(OBJS) libskit.a libblas.a
#	$(F90) -o ptetra $(FLAGS) $(OBJS) $(SPARSKIT) $(BLAS)

ptetra.o:       ptetra.f90
	$(F90) -c $(FLAGS) ptetra.f90

clean:
	rm -f ptetra $(OBJS) *.mod *__*

