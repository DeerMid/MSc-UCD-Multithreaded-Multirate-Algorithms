#include<iostream>
#include<vector>
#include "SRsolvers.h"
#include "systems.h"
#include "MRsolvers.h"
#include<stdlib.h>

using namespace std;



int main(int argc, char **argv) {

	//default tolerance values
	const float relTol = 0.1;
	const float absTol = 0.01;
	cout.precision(4);

	//check user input
	if(argc != 5){
		
		cout << "Please input the following arguments when calling the program: Start Time, End Time, Step Size, Dimension" << endl;
		return 1;
	}
	else{

		//initial conditions
		float tStart = atof(argv[1]);
		vector<float> t = { tStart };
		vector<float> yRand;
		float dim = atoi(argv[4]);
		for(int i = 0; i < dim; i++) yRand.push_back(static_cast <float> (rand()) / static_cast <float> (RAND_MAX/5.0));
		float h = atof(argv[3]);
		float T = atof(argv[2]);
		
		vector<float> y = yRand;
		
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
		for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim + 2) << endl;
		}

		cout << "#############################" << endl << endl;
		cout << "RK12 Multirate Method Results" << endl;

		t = { tStart };

		y = yRand;

		internalInfo inter;
		solRK12MR(t, y, h, T, &linearScalarODE, relTol, absTol, dim, inter);

		cout << "\nt\ty" << endl;
		cout << "---------------" << endl;
		for (int i = 0; i < t.size(); i++) {
			cout << t.at(i) << "\t" << y.at(i * dim + 2) << endl;
		}

		return 0;
	}

}
