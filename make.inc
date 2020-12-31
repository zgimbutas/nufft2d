ifeq ($(HOST),linux)

OBJSUF=o
MODSUF=mod
FC=f77 -c 
FFLAGS=-fast
FLINK=f77 -o $(PROJECT)

endif

ifeq ($(HOST),macosx)

OBJSUF=o
MODSUF=mod
FC=gfortran -c 
FFLAGS=-O3
FLINK=gfortran -o $(PROJECT)

endif


ifeq ($(HOST),linux-openmp)

OBJSUF=o
MODSUF=mod
FC=f77 -c
FFLAGS=-fast -fopenmp
FLINK=f77 -fopenmp -o $(PROJECT)

endif

ifeq ($(HOST),macosx-openmp)

OBJSUF=o
MODSUF=mod
FC=gfortran -c -fopenmp
FFLAGS=-O3
FLINK=gfortran -fopenmp -o $(PROJECT)

endif


