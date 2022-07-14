#include<iostream>
#include<vector>
#include "SRsolvers.h"
#include "systems.h"

using namespace std;



int main() {

	//initial conditions
	vector<float> t = { 0.0 };
	vector<float> y = { 1.0 };

	solForwardEulerSR(t, y, 0.78, 5, &linearScalerODE); //compute forward Euler

	//print results
	cout << "Forward Euler Method Results" << endl;

	cout << "\nt\ty" << endl;
	cout << "---------------" << endl;
	for (int i = 0; i < t.size(); i++) {
		cout << t.at(i) << "\t" << y.at(i) << endl;
	}

	cout << "#############################" << endl << endl;
	cout << "Heun Method Results" << endl;

	t = { 0.0 };
	y = { 1.0 };

	solHeunSR(t, y, 0.78, 5, &linearScalerODE); //compute for Heun's

	cout << "\nt\ty" << endl;
	cout << "---------------" << endl;
	for (int i = 0; i < t.size(); i++) {
		cout << t.at(i) << "\t" << y.at(i) << endl;
	}

	cout << "#############################" << endl << endl;
	cout << "RK12 Method Results" << endl;

	t = { 0.0 };
	y = { 1.0 };

	solRK12SR(t, y, 0.78, 5, &linearScalerODE, 0.1, 0.01); //compute for RK12

	cout << "\nt\ty" << endl;
	cout << "---------------" << endl;
	for (int i = 0; i < t.size(); i++) {
		cout << t.at(i) << "\t" << y.at(i) << endl;
	}



	return 0;
}