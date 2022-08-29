//this file contains implementations of the singlerate **single-dimension** solvers in C++

#include "SRsolvers.h"
#include<iostream>
#include<vector>
#include<math.h>

using namespace std;

void solForwardEulerSR(vector<float>& t, vector<float>& y, float h, float T, float (*f)(float, float), int dim) {

	//the basic forward Euler implmentation
	// vector<float> t represents the time vector, containing each time step over the interval
	// vector <float> y represents the computed values at the corresponding time step
	// int dim accounts for striding of the y vector

	float t0 = t.back();
	int count = 0;

	while (t0 < T) {
		for (int extend = 0; extend < dim; extend++) y.push_back(0.0); //initialise the next "row" of the flattened 2D vector
		for (int i = 0; i < dim; i++) y[i + (dim * (count + 1))] = y[i + (dim * count)] + h * f(t0, y[i + (dim * count)]); //forward Euler
		t0 = t0 + h;
		t.push_back(t0);

		if (t0 + h > T) { //this aligns the endpoint so we don't exceed the required interval
			h = T - t0;
		}
		count++; //iterate the count variable to force computations to be conducted on the next row for the following loop
	}
}


void solHeunSR(vector<float>& t, vector<float>& y, float h, float T, float (*f)(float, float), int dim) {

	//implementation of Heun's method
	// vector<float> t represents the time vector, containing each time step over the interval
	// vector <float> y represents the computed values at the corresponding time step
	// int dim accounts for striding of the y vector

	float t0 = t.back();
	int count = 0;

	while (t0 < T) {
		for (int extend = 0; extend < dim; extend++) y.push_back(0.0); //initialise the next "row" of the flattened 2D vector
		for (int i = 0; i < dim; i++) {
			float k1 = h * f(t0, y[i + (dim * count)]); //compute both components of Heun's
			float k2 = h * f(t0 + h, y[i + (dim * count)] + k1);

			y[i + (dim * (count + 1))] = y[i + (dim * count)] + 0.5 * (k1 + k2); //computes and assigns the Heun computation to the next "vertical entry" 
		}
		t.push_back(t0 + h);
		t0 = t0 + h;

		if (t0 + h > T) {
			h = T - t0;
		}
		count++; //iterate the count variable to force computations to be conducted on the next row for the following loop
	}
}

void solRK12SR(vector<double>& t, vector<double>& y, double h, double T, void (*f)(double, vector<double>&, int, int, vector<double>&, vector<bool>&), double relTol, double absTol, int dim) {

    //implementation of the RK12 method with adaptive time step
    // vector<float> t represents the time vector, containing each time step over the interval
    // vector <float> y represents the computed values at the corresponding time step
    // int dim accounts for striding of the y vector

    //these constants determine the time step refinement and error bounds
    const double nu = 0.7;
    const double hLo = 0.1;
    const double hHi = 5;
    double fac = 1.0; //preinitialise

    double t0 = t.back();
    int count = 0;

    vector<double> k1 (dim, 0.0);
    vector<double> k2 (dim, 0.0);
    vector<double> nErrVec (dim, 0.0);
    vector<bool> dumRef (dim, true);

    while (t0 < T) {
	double nErrMax = 0.0;
        for (int i = 0; i < dim; i++) {
		vector<bool> dumRef (dim, true);
                f(t0, y, dim, count, k1, dumRef); //compute both components of Heun's
                for(int i = 0; i < dim; i++) k1[i] = k1[i] * h;
		vector<double> k2temp (dim, 0.0);
                for(int i = 0; i < dim; i++) k2temp[i] = y[i + (dim * count)] + k1[i];
                f(t0 + h, k2temp, dim, 0, k2, dumRef);
                for(int i = 0; i < dim; i++) k2[i] = h * k2[i];
	}
	for(int i = 0; i < dim; i++){
        	double err = 0.5 * (k2[i] - k1[i]);
                err = fabs(err); //compute error comparing forward and Heun

                double tolC = relTol * fabs(y[i + (dim * count)]) + absTol;
                double nErr = err / tolC;
                nErrVec[i] = nErr; //stores that node's normalised error
	}
	for(int i = 0; i < dim; i++) if(nErrMax < nErrVec[i]) nErrMax = nErrVec[i]; 
	fac = nu * pow(nErrMax, -0.5);

        if (nErrMax > 1.0) { //if the error is too large we rescale the timestep
                fac = max(hLo, fac);
                h = fac * h;
        }
        else { //otherwise continue with appending
		for (int extend = 0; extend < dim; extend++) y.push_back(0.0);

                for(int i = 0; i < dim; i++) y[i + (dim * (count + 1))] = y[i + (dim * count)] + 0.5 * (k1[i] + k2[i]); //computes and assigns the computation to the next "vertical entry"
        
        	t0 = t0 + h;
        	t.push_back(t0);

        	fac = min(hHi, fac); //reset fac
        	h = fac * h;

        	if (t0 + h > T) {
            		h = T - t0;
        	}
		count++; //iterate the count variable to force computations to be conducted on the next row for the following loop

	}	
    }
}

