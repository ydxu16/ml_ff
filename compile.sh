#!/bin/bash

# MKL compile/link flags, see mkl link advisor
MKLROOT="/home/kuangy/opt/intel/mkl/"
MKL_CFLAGS="-DMKL_ILP64 -m64 -I${MKLROOT}/include"
MKL_LDFLAGS="-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"

# user defined compile/link flags
CFLAGS="-I/home/kuangy/data2/data/project_fhc_ml/playground/cpp/include"
LDFLAGS="-L/home/kuangy/data2/data/project_fhc_ml/playground/cpp/lib"
# user defined libs
LIBS=""


CFLAGS="$CFLAGS $MKL_CFLAGS"
LDFLAGS="$LDFLAGS $MKL_LDFLAGS"

cd lib
rm *so
sources=`ls *cpp`
for src in $sources; do
    lib=${src/.cpp/}
    LIBS="$LIBS -l$lib"
    g++ -fPIC $lib.cpp -o lib$lib.so -shared $CFLAGS $LDFLAGS
done
cd ..
LDFLAGS="$LDFLAGS $LIBS"

g++ test.cpp -o main $CFLAGS $LDFLAGS
