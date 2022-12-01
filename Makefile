CFLAGS=-O3 -g -W -Wall -Wextra -ffast-math
CXXFLAGS=$(CFLAGS) -Wno-unused

simulator: simulator.o cross-correlation.o terminal-canvas.o
	g++ -std=c++11 -O3 -W -Wall -Wextra $^ -o $@ 

simulator.o: cross-correlation.h
cross-correlation.o: cross-correlation.h

clean:
	rm -f *.o simulator
