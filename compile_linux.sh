#!/bin/sh
#
#  compile_linux.sh - build the Octave (.mex) and MATLAB (.mexa64) MEX
#  files of nufft2d on Linux x86-64 from nufft2d.c (mwrap output of
#  nufft2d.mw) and the F90 Fortran sources.
#
#  Octave needs mkoctfile (Debian/Ubuntu package octave-dev, Fedora
#  octave-devel). The MATLAB file is built when MATLAB is set to the
#  MATLAB root: its bin/mex compiles the gateway and links the gfortran
#  objects with the system libgfortran.
#
#      sh compile_linux.sh
#      MATLAB=/usr/local/MATLAB/R2024b sh compile_linux.sh
#
#  FC and MKOCTFILE override the compilers.
#
set -e
P=nufft2d
SRC="nufft2df90.f dirft2d.f dfft.f next235.f"
FC=${FC:-gfortran}
MKOCTFILE=${MKOCTFILE:-mkoctfile}
B=$(mktemp -d)
cp $SRC $P.c $B
cd $B
F77=$FC FFLAGS="-O3 -std=legacy -fPIC -w" $MKOCTFILE --mex -o $P.mex $P.c $SRC $($MKOCTFILE -p FLIBS)
if [ -n "$MATLAB" ]; then
  rm -f *.o
  for f in $SRC; do $FC -O3 -std=legacy -fPIC -w -c $f; done
  $MATLAB/bin/mex -largeArrayDims $P.c *.o -lgfortran -lm
fi
cd - > /dev/null
cp $B/$P.mex* .
rm -rf $B
