all:
	g++ SRsolvers.cpp MRsolvers.cpp systems.cpp main.cpp -o solvers -O0 -g -std=c++11 

clean:
	rm -f solvers
