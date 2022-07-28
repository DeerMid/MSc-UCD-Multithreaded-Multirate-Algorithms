all:
	g++ SRsolvers.cpp MRsolvers.cpp systems.cpp main.cpp -o srsolvers -O0 -g -std=c++11 

clean:
	rm -f srsolvers
