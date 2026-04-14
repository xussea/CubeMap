FC=gfortran
#FFLAGS=-O2 -fopenmp -m64 -I${MKLROOT}/include -fcheck=all -fbacktrace -Wall -Wextra -pedantic
#LIBFLAGS=-L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl
FFLAGS=-O2 -fdefault-integer-8 -fopenmp -m64 -I${MKLROOT}/include -fcheck=all -fbacktrace -Wall -Wextra -pedantic
LIBFLAGS=-Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl
OUT=cubemap

all: cubemap clean

cubemap:
	${FC} ${FFLAGS} CubeMap.f90  ${LIBFLAGS} -o ${OUT}
		mkdir -p ../bin
			mv ${OUT} ../bin
clean:
	rm -f *.o *.mod
