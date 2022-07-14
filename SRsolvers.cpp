//this file contains implementations of the singlerate **single-dimension** solvers in C++


#include<iostream>
#include<vector>

using namespace std;

void solForwardEulerSR(vector<float>& t, vector<float>& y, float h, float T, float (*f)(float, float)) {
	
	//the basic forward Euler implmentation
	// vector<float> t represents the time vector, containing each time step over the interval
	// vector <float> y represents the computed values at the corresponding time step

	float y0 = y.back(); //store initial values
	float t0 = t.back();

	while (t0 < T) {
		y0 = y0 + h * f(t0, y0); //forward Euler
		y.push_back(y0); //update
		t0 = t0 + h;
		t.push_back(t0);

		if (t0 + h > T) { //this aligns the endpoint so we don't exceed the required interval
			h = T - t0;
		}
	}
}

void solHeunSR(vector<float>& t, vector<float>& y, float h, float T, float (*f)(float, float)) {

	//implementation of Heun's method
	// vector<float> t represents the time vector, containing each time step over the interval
	// vector <float> y represents the computed values at the corresponding time step

	float y0 = y.back();
	float t0 = t.back();

	while (t0 < T) {
		float k1 = h * f(t0, y0); //compute both components of Heun's
		float k2 = h * f(t0 + h, y0 + k1);

		y0 = y0 + 0.5 * (k1 + k2); //y0 represents the next value at t0+h

		y.push_back(y0); //update vectors with latest time step
		t0 = t0 + h;
		t.push_back(t0);

		if (t0 + h > T) {
			h = T - t0;
		}
	}
}

void solRK12SR(vector<float>& t, vector<float>& y, float h, float T, float (*f)(float, float), float relTol, float absTol) {

	//implementation of the RK12 method with adaptive time step
	// vector<float> t represents the time vector, containing each time step over the interval
	// vector <float> y represents the computed values at the corresponding time step

	//these constants determine the time step refinement and error bounds
	const float nu = 0.7;
	const float hLo = 0.1;
	const float hHi = 5;

	float y0 = y.back();
	float t0 = t.back();

	while (t0 < T) {
		float k1 = h * f(t0, y0);
		float k2 = h * f(t0 + h, y0 + k1);

		float err = 0.5 * (k1 - k2);
		err = abs(err); //compute error comparing forward and Heun's

		float tolC = relTol * abs(y0) + absTol;
		float nErr = err / tolC;

		float fac = nu * pow(nErr, -0.5);

		if (nErr > 1.0) { //if the error is too large we rescale the timestep
			fac = max(hLo, fac);
			h = fac * h;
		}
		else { //otherwise continue
			y0 = y0 + 0.5 * (k1 + k2); //y0 represents the next value at t0+h

			y.push_back(y0); //update vectors with latest time step
			t0 = t0 + h;
			t.push_back(t0);

			fac = min(hHi, fac);
			h = fac * h;

			if (t0 + h > T) {
				h = T - t0;
			}
		}
		}
	}
