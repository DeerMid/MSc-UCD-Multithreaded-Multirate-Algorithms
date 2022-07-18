#include<iostream>
#include<vector>
#include "SRsolvers.h"
#include "systems.h"

using namespace std;



int main() {

	//initial conditions
	vector<float> t = { 0.0 };
	vector<float> y = { 1.0 , 2.0, 3.0};
	int dim = 3;

	solForwardEulerSR(t, y, 0.78, 5, &linearScalerODE, dim); //compute forward Euler

	//print results
	cout << "Forward Euler Method Results" << endl;

	cout << "\nt\ty" << endl;
	cout << "---------------" << endl;
	for (int i = 0; i < t.size(); i++) {
		cout << t.at(i) << "\t" << y.at(i * dim) << endl; //output the Euler method for the first node point
	}

	cout << "#############################" << endl << endl;
	cout << "Heun Method Results" << endl;

	t = { 0.0 };
	y = { 1.0 , 2.0, 3.0 };

	solHeunSR(t, y, 0.78, 5, &linearScalerODE, dim); //compute for Heun's method for the second node point

	cout << "\nt\ty" << endl;
	cout << "---------------" << endl;
	for (int i = 0; i < t.size(); i++) {
		cout << t.at(i) << "\t" << y.at(i * dim + 1) << endl;
	}

	cout << "#############################" << endl << endl;
	cout << "RK12 Method Results" << endl;

	t = { 0.0 };
	y = { 1.0 , 2.0, 3.0 };

	solRK12SR(t, y, 0.78, 5, &linearScalerODE, 0.1, 0.01, dim); //compute for RK12 method for the third node point

	cout << "\nt\ty" << endl;
	cout << "---------------" << endl;
	for (int i = 0; i < t.size(); i++) {
		cout << t.at(i) << "\t" << y.at(i * dim + 2) << endl;
	}



	return 0;
}