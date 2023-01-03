CXX?=g++
NVCC?=nvcc
CFLAGS=-O3 -W -Wall -Wextra
EXTRA_LDFLAGS?=
EXTRA_CCFLAGS?=
CXXFLAGS=$(CFLAGS) -Wno-unused -std=c++20 $(EXTRA_LDFLAGS) $(EXTRA_CCFLAGS)
#CUDA_OBJECTS=

all: simulator cuda-fft-example

simulator: simulator.o cross-correlation.o terminal-canvas.o
	$(CXX) $(CXXFLAGS) $^ -o $@ -lfftw3

simulator.o: cross-correlation.h
cross-correlation.o: cross-correlation.h

%.o: %.cu
	$(NVCC) -c $< -o $@

# cuda-backend.o: $(CUDA_OBJECTS)
# 	$(NVCC) –dlink $^ –o $@

cuda-fft-example: cuda-fft-example.o
	$(CXX) $(CXXFLAGS) cuda-fft-example.o -o $@ -lcudart

clean:
	rm -f *.o simulator cuda-fft-example
