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
	const float deltaRK = 5e-1;
	float fac = 1.0; 
	float nErrMax = 1.0; //preinitialise
	float nref_nErrMax = 0.0;

	vector<float> nErrVec; //vector to store a time step's normalised error	

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
				nErrVec[i] = nErr; //stores that node's normalised error
			}

			//check for the maximum error
			for(int i = 0; i < dim; i++){
				if(nErrMax <= nErrVec[i]) nErrMax = nErrVec[i];
			}

			vector<bool> ref (dim, false); //create a default refinement vector where no node needs refining
			
			//now we check for refinement
			if(nErrMax > 1.0){ //if the worst error is outside tolerance then some refinement is necessary
				for(int i = 0; i < dim; i++){
					if(nErrVec[i] > (deltaRK*nErrMax)) ref[i] = false; //identifies nodes which need refinement
				}
			}

			//we now need to identify the max error of the NOT to be refined set
			for(int i = 0; i < dim; i++){
				if(ref[i] == false){
					if(nref_nErrMax < nErrVec[i]) nref_nErrMax = nErrVec[i];
				}
			}
			
			//this step determines how "large" the refinement will be, this can probably be improved much further in terms of efficiency
			int check = 0;
			for(int i = 0; i < dim; i++) check += ref[i];

			if(check == dim){ //all components need to be refined
				fac = nu * pow(nErrMax, -0.5);
				fac = max(hLo, fac);
				h = fac * h;
			}
			else if(nref_nErrMax > 1.0){
				fac = nu * pow(nref_nErrMax, -0.5);
				fac = max(hLo, fac);
				h = fac * h;
			}
			else{
				if(check > 0){ //this checks if there are components that need refinement
					float Tin = t0 + h;
					fac = nu * pow(nErrMax, -0.5);
					fac = max(hLo, fac);
					float hin = fac * h;

					//interpolation variables
					vector<float> yInter;
					//the following line equates to yInter = y + 0.5 * (k1 + k2) from above
					//maybe it's more efficient to store the k1 and k2 variables as vectors and recall here
					for(int i = 0; i < dim; i++) yInter[i] = y[i + (dim * inter.count)] + 0.5 * ((h * f(t0 + h, y[i + (dim * inter.count)]) + (h * f(t0, y[i + (dim * inter.count)])));
					vector<float> y1nref;
					for(int i = 0; i < dim; i++) y1nref[i] = yInter[i] * (int) nref[i];
				}
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
