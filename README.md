# CubeMap approach
## Instalation
You run the CubeMap code at your own risk and you are responsible of checking that the obtained data makes sense.

**Compile the code**
Before compiling, make sure you have installed:

Intel MKL 2025
GNU Fortran (gfortran) 15

Also ensure that the environment variable MKLROOT is correctly defined.

  1. gfortran -O2 -fdefault-integer-8 -fopenmp -m64 -I${MKLROOT}/include -fcheck=all -fbacktrace -Wall -Wextra -pedantic CubeMap.f90 -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl -o cubemap
  2. Using the Makefile of the src/ folder:
```bash
make
```

**Run the code**
You need:
 1. A cube file of the molecular orbital of interest (generated with any quantum chemistry software)
 2. The .xyz files of the two molecules of interest to compute the coupling

It is runned as:

````bash
cubemap --c cubefile --m xyz1 xyz2 
````

**Excluding atoms**
To exclude atoms from the CubeMap calculation, create an Exclude file containing:
 1. First line: number of atoms to exclude
 2. Second line: index of those atoms

For example:
````
10
1 3 4 6 7 12 13 16 20 21
````


It is runned as:
````bash
cubemap --c cubefile --m xyz1 xyz2 --ex Exclude
````

**Multiple Cubes Files**
This code also supports using multiple reference cube files, with optional atom exclusion:

````bash
cubemap --mc cubefile1 cubefile2 xyzref1 xyzref2 box --m xyz1 xyz2 [--ex Excludefile]
````

**Help**
You can display all available options of cubemap by running:

````bash
cubemap --h [--help]
````

Disclaimer: You run the CubeMap code at your own risk and you are responsible of checking that the obtained data makes sense.
