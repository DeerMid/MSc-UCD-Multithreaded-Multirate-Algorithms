all:
	g++ double_MRsolvers.cpp double_HEode.cpp SRsolvers.cpp double_main2.cpp -o solversNew_double -O0 -g -std=c++11
	g++ double_MRsolvers.cpp OMPdouble_HEode.cpp SRsolvers.cpp OMPdouble_main2.cpp -o OMPsolversNew_double -O0 -g -std=c++11 -fopenmp

clean:
	rm -f solversNew_double OMPsolversNew_double
