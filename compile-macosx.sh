
PATH=/Applications/MATLAB_R2011b.app/bin/:$PATH

make -f makefile.mwrap -j4 TARGET=matlab-maci64 clean
make -f makefile.mwrap -j4 TARGET=matlab-maci64  

make -f makefile.mwrap clean distclean

