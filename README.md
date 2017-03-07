# Nekbone_kernel with LIBXSMM/SIMD

All the compilers/flags and defination in test/makenek

### General information

```
# Fortran compiler
F77="mpiifort"

# C compiler
CC="mpiicc"
```

### Using the in-hourse SIMD implementation

```
uncomment out the variable "IFSIMD"
# Enable SIMD (default false)
#IFSIMD="true"
```

### Using the libxsmm

uncomment out the variable "IFXSMM"

```
# Enable XSMM (default false)
#IFXSMM="true"
```
In the case, include paths and linked libraried should be set

```
USR_INCDIR="..."
USR_INCDIR="..."
```

### the polynomial order and the number of elements per core can
be set in fiel SIZE

### Compiling and running the code

```
$ cd test
$ ./makenek clean
$ ./makenek test
$ mpirun -n  ./nekkernel
```




