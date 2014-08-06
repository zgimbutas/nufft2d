
make -f makefile.mwrap -j4 TARGET=matlab-windows-w32 clean
make -f makefile.mwrap -j4 TARGET=matlab-windows-w32 

make -f makefile.mwrap -j4 TARGET=matlab-windows-w64 clean
make -f makefile.mwrap -j4 TARGET=matlab-windows-w64 

make -f makefile.mwrap -j4 TARGET=matlab-linux-a64 clean
make -f makefile.mwrap -j4 TARGET=matlab-linux-a64  

make -f makefile.mwrap -j4 TARGET=octave-linux clean
make -f makefile.mwrap -j4 TARGET=octave-linux 

make -f makefile.mwrap clean distclean

