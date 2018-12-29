
CC=g++

# USE MKL or eigen/omp
#MKLROOT=/home/kuangy/opt/intel/mkl/
#MKL_CFLAGS=-DMKL_ILP64 -m64 -I$(MKLROOT)/include
#MKL_LDFLAGS=-L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
MKLROOT=
MKL_CFLAGS=
MKL_LDFLAGS=
USEOMP= -fopenmp

# user defined compile/link flags
COPTS= -g -O3 -Wall -std=c++11 -mavx2 -I/home/kuangy/data2/data/project_fhc_ml/playground/cpp/include $(USEOMP)
LDOPTS= -L/home/kuangy/data2/data/project_fhc_ml/playground/cpp/lib
PY_INCLUDE= -I/usr/include/python3.6m -I/usr/local/include/python3.6 -I/home/kuangy/.local/include/python3.6m

CFLAGS=$(COPTS) $(MKL_CFLAGS) $(PY_INCLUDE)
LDFLAGS=$(LDOPTS) $(MKL_LDFLAGS)

# user defined libs
LIBS= lib/libmylib.so
OBJS= lib/file_io.o \
	  lib/spatial.o \
	  lib/neighbour_list.o
PYBDERS = spatial.so neighbour_list.so
LIBLINKS= -lmylib -lpython3.6m
HEADERS= include/mylibs.h \
		 include/constants.h \
		 include/spatial.h \
		 include/file_io.h \
		 include/neighbour_list.h
TARGET = main

default: all

all: $(TARGET) $(OBJS) $(LIBS)

python: $(PYBDERS)

$(TARGET): test.cpp $(LIBS) $(HEADERS) $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) test.cpp $(LIBLINKS)

lib/file_io.o: lib/file_io.cpp include/file_io.h
	$(CC) -fPIC lib/file_io.cpp -c -o lib/file_io.o $(CFLAGS) $(LDFLAGS)

lib/spatial.o: lib/spatial.cpp include/spatial.h
	$(CC) -fPIC lib/spatial.cpp -c -o lib/spatial.o $(CFLAGS) $(LDFLAGS)

lib/neighbour_list.o: lib/neighbour_list.cpp include/neighbour_list.h
	$(CC) -fPIC lib/neighbour_list.cpp -c -o lib/neighbour_list.o $(CFLAGS) $(LDFLAGS)

lib/libmylib.so: $(OBJS)
	$(CC) -fPIC $(OBJS) -shared -o lib/libmylib.so $(CFLAGS) $(LDFLAGS)

clean: 
	$(RM) lib/*so *so lib/*o main

# python3 binders
spatial.so: lib/spatial.cpp
	g++ -O3 -Wall -shared -fPIC lib/spatial.cpp -o spatial.so $(CFLAGS) $(LDFLAGS)

neighbour_list.so: lib/neighbour_list.cpp lib/spatial.cpp
	g++ -O3 -Wall -shared -fPIC lib/spatial.cpp lib/neighbour_list.cpp -o neighbour_list.so $(CFLAGS) $(LDFLAGS)

