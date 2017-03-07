# nekbone_kernel

All the compilers/flags and defination in test/makenek

1) General information

# Fortran compiler
F77="mpiifort"

# C compiler
CC="mpiicc"

2) Using the in-hourse SIMD implementation

uncomment out the variable "IFSIMD"
# Enable SIMD (default false)
#IFSIMD="true"


3) Using the libxsmm

uncomment out the variable "IFXSMM"
# Enable XSMM (default false)
#IFXSMM="true"

In the case, include paths and linked libraried should be set

USR_INCDIR="..."
USR_INCDIR="..."

4) the polynomial order and the number of elements per core can
be set in fiel SIZE

5) Compiling and running the code

$ cd test
$ ./makenek clean
$ ./makenek test
$ mpirun -n  ./nekkernel





