#!/bin/sh
#
#  compile_macos.sh - build the Octave (.mex) and MATLAB (.mexmaca64) MEX
#  files of nufft2d on an Apple Silicon Mac from nufft2d.c (mwrap output of
#  nufft2d.mw) and the F90 Fortran sources.
#
#  Octave: Homebrew octave, built with Homebrew gcc; use its gfortran so
#  that one libgfortran is loaded. MATLAB: mex (Xcode clang) on Fortran
#  objects from the same gfortran. libgfortran is linked through
#  Homebrew's stable opt/gcc/.../current path.
#
#      sh compile_macos.sh            (FC, MATLAB to override)
#
set -e
P=nufft2d
SRC="nufft2df90.f dirft2d.f dfft.f next235.f"
FC=${FC:-/opt/homebrew/opt/gcc/bin/gfortran}
MATLAB=${MATLAB:-/Applications/MATLAB_R2024b.app}
FL=/opt/homebrew/opt/gcc/lib/gcc/current
B=$(mktemp -d)
cp $SRC $P.c $B
cd $B
F77=$FC FFLAGS="-O3 -std=legacy -fPIC -w" mkoctfile --mex -o $P.mex $P.c $SRC $(mkoctfile -p FLIBS)
rm -f *.o
for f in $SRC; do $FC -O3 -std=legacy -fPIC -w -c $f; done
$MATLAB/bin/mex -largeArrayDims $P.c *.o -L$FL -lgfortran -lquadmath LDFLAGS='$LDFLAGS -Wl,-rpath,'$FL
cd - > /dev/null
cp $B/$P.mex $B/$P.mexmaca64 .
rm -rf $B
