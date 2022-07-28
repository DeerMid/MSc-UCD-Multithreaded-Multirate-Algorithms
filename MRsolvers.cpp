#include "MRsolvers.h"
#include<vector>
#include<iostream>
#include<math.h>

using namespace std;

void solRK12MR(vector<float>& t, vector<float>& y, float h, float T, float (*f)(float, float), float relTol, float absTol, int dim, internalInfo inter) {

	//these constants determine the time step refinement and error bounds
	const float nu = 0.7;
	const float hLo = 0.1;
	const float hHi = 5;
	float fac = 1.0; 
	float nErrMax = 0.0; //preinitialise
	
	float t0 = t.back();

	while (t0 < T) {
		if (inter.internalCheck == false) { //we are in the outer loop
			for (int extend = 0; extend < dim; extend++) y.push_back(0.0);
			for (int i = 0; i < dim; i++) {
				float k1 = h * f(t0, y[i + (dim * inter.count)]); //compute both components of Heun's
				float k2 = h * f(t0 + h, y[i + (dim * inter.count)] + k1);

				float err = 0.5 * (k1 - k2);
				err = fabs(err); //compute error comparing forward and Heun's

				float tolC = relTol * fabs(y[i + (dim * inter.count)]) + absTol;
				float nErr = err / tolC;
				if (nErrMax < nErr) nErrMax = nErr; //stores the largest error value

			}
			//debugging tool, force closes the while loop and prints a message 
			cout << "Made it this far!" << endl;
			t0 = T + 100;
		}
		else {
			cout << "This is the internal refinement stage, come back later!" << endl;
		}
	}
}
