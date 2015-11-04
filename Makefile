all: spsim
clean:
	rm spsim *.o

#                       NOTE:
# MPI = ON     compile with MPI support
# MPI = OFF    no MPI support

MPI = OFF


LIBCONFIG = libconfig-0.9
CFLAGS =   -I$(LIBCONFIG)

FLAGS = -O3

ifeq ($(MPI),ON)
CC = mpicc -Wall
CXX = mpic++ -Wall
CFLAGS += -DMPI
else
CC = gcc -Wall $(FLAGS)
CXX = g++ -Wall $(FLAGS) -std=c++11
CU = /usr/local/cuda/bin/nvcc -arch=sm_20 -lineinfo -I /usr/local/cuda/include/ $(FLAGS)
endif

LDFLAGS = $(FLAGS)
LOADLIBES =  -lm -L/usr/local/cuda/lib64 -lcuda -lcudart

cuda_compute_pattern.o: cuda_compute_pattern.cu cuda_compute_pattern.h
	$(CU) -c cuda_compute_pattern.cu -o $@

LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)


spsim: spsim.o config.o diffraction.o molecule.o io.o mpi.o noise.o amplification.o random_rotation.o random_position.o cuda_compute_pattern.o $(LIBCONFIG)/.libs/libconfig.a


