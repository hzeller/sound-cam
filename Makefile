CFLAGS=-O3 -g -W -Wall -Wextra
CXXFLAGS=-O3 -W -Wall -Wextra -Wno-unused -ffast-math -DNDEBUG

simulator: simulator.o cross-correlation.o terminal-canvas.o
	g++ -std=c++11 -O3 -W -Wall -Wextra $^ -o $@ -lalglib3

simulator.o: cross-correlation.h
cross-correlation.o: cross-correlation.h

clean:
	rm -f *.o simulator
