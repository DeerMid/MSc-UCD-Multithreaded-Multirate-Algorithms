#include<iostream>
#include<vector>
#include "double_HEode.h"
#include "double_MRsolvers.h"
#include<stdlib.h>
#include<math.h>

using namespace std;
#define _USE_MATH_DEFINES


int main(int argc, char **argv) {

	//default tolerance values
	const double relTol = 0.0;
	const double absTol = 1e-2;
//	cout.precision(4);

	//check user input
	if(argc != 5){
		
		cout << "Please input the following arguments when calling the program: Start Time, End Time, Step Size, Dimension" << endl;
		return 1;
	}
	else{

		//initial conditions
		double tStart = atof(argv[1]);
		vector<double> t = { tStart };
		//vector<double> yRand;
		int dim = atoi(argv[4]);
		vector<double> y (dim, 0.0);
		//for(int i = 0; i < dim; i++) yRand.push_back(static_cast <double> (rand()) / static_cast <float> (RAND_MAX/5.0));
		double h = atof(argv[3]);
		double T = atof(argv[2]);
		//cout << dim << endl;	
	//	vector<double> y = yRand;
	/*	
		solForwardEulerSR(t, y, h, T, &linearScalarODE, dim); //compute forward Euler

		//print results
		cout << "Forward Euler Method Results" << endl;

		cout << "\nt\ty" << endl;
		cout << "---------------" << endl;
		for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim) << endl; //output the Euler method for the first node point
		}

		cout << "#############################" << endl << endl;
		cout << "Heun Method Results" << endl;

		t = { tStart };

		y = yRand;

		solHeunSR(t, y, h, T, &linearScalarODE, dim); //compute for Heun's method for the second node point

		cout << "\nt\ty" << endl;
		cout << "---------------" << endl;
		for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim + 1) << endl;
		}
	
		cout << "#############################" << endl << endl;
		cout << "RK12 Method Results" << endl;

		t = { tStart };

		y = yRand;

		solRK12SR(t, y, h, T, &linearScalarODE, relTol, absTol, dim); //compute for RK12 method for the third node point

		cout << "\nt\ty" << endl;
		cout << "---------------" << endl;
		/*for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim + 2) << endl;
		}
		*/
		//cout << "The size is " << t.size() << endl;
		cout << "#############################" << endl << endl;
		cout << "RK12 Multirate Method Results" << endl;

		t = { tStart };

		//y = yRand;
		for(int i = 0; i < dim; i++) y[i] = (1.0*i)/(dim - 1);
		for(int i = 0; i < dim; i++) cout << " y[" << i << "] = " << y[i] << " ," << endl;  
		for(int i = 0; i < dim; i++) y[i] = (double) sin(2 * M_PI * y[i]);
		for(int i = 0; i < dim; i++) cout << " y[" << i << "] = " << y[i] << " ," << endl;  
		
		internalInfo inter;
		int rowCount = 0;
		solRK12MR(t, y, h, T, &HEode, relTol, absTol, dim, inter, rowCount);

		cout << "\nt\ty[0]" << endl;
		cout << "---------------" << endl;
		for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim) << endl;
		}

		
		cout << "\nt\ty[1]" << endl;
		cout << "---------------" << endl;
		for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim + 1) << endl;
		}


		cout << "\nt\ty[2]" << endl;
		cout << "---------------" << endl;
		for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim + 2) << endl;
		}
		/*
		cout << "\nt\ty[37]" << endl;
		cout << "---------------" << endl;
		for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim + 37) << endl;
		}
		*/
		cout << "final t-1 is " << t.at(t.size()-2);
		cout << ", y-1 is " << y.at(t.size()-2) << endl;

		cout << "final t is " << t.at(t.size()-1);
		cout << ", y is " << y.at(t.size()-1) << endl;	
		cout << "The size is " << t.size() << endl;

		cout << "test: 3.0 *t.size()/y.size() = " << 3.0 * t.size()/y.size() << endl;
		//for(int i=0; i<y.size(); i++) cout << "i: " << i << ", y:" <<  y.at(i) <<endl;
		//cout << endl << "Row count is: " << rowCount << endl;
		cout << "the size of t is " << t.size() << endl;	
		return 0;
	}

}
