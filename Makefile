CFLAGS=-O3 -g -W -Wall -Wextra
CXXFLAGS=-O3 -W -Wall -Wextra -Wno-unused

simulator: simulator.o cross-correlation.o terminal-canvas.o
	g++ -std=c++11 -O3 -W -Wall -Wextra $^ -o $@

clean:
	rm -f *.o simulator
