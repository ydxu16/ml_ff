#!/bin/bash
MKLLIB=/home/kuangy/opt/intel/mkl/lib/intel64/
export LD_PRELOAD=${MKLLIB}/libmkl_def.so:${MKLLIB}/libmkl_avx2.so:${MKLLIB}/libmkl_core.so:${MKLLIB}/libmkl_intel_lp64.so:${MKLLIB}/libmkl_gnu_thread.so
export LD_PRELOAD=$LD_PRELOAD:/usr/lib/gcc/x86_64-linux-gnu/7/libgomp.so
export OMP_NUM_THREADS=12

python3 -m cProfile -s cumtime ./pytest.py > log1
