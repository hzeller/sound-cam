simulator: simulator.o terminal-canvas.o
	g++ -std=c++11 -O3 -Wall -Wextra $^ -o $@
