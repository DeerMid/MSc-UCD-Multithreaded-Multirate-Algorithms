#include "double_MRsolvers.h"
#include<vector>
#include<iostream>
#include<math.h>

using namespace std;

void solRK12MR(vector<double>& t, vector<double>& y, double h, double T, void (*f)(double, vector<double>&, int, int, vector<double>&, vector<bool>&), double relTol, double absTol, int dim, internalInfo &inter, int &rowCount) {

	//these constants determine the time step refinement and error bounds
	const double nu = 0.7;
	const double hLo = 0.1;
	const double hHi = 5;
	const double deltaRK = 5e-1;
	const double nuG = 1.0;
	const double hHiG = 25.0;
	double fac = 1.0; 

	double nErrMax = 0.0; //preinitialise error maximums
	double nref_nErrMax = 0.0;

	//create vectors to contain the components of Heun's method
	vector<double> k1 (dim, 0.0);
	vector<double> k2 (dim, 0.0);

	vector<double> nErrVec (dim, 0.0); //vector to store a time step's normalised error	

	double t0 = t.back();

	while (t0 < T) {
		if (inter.internalCheck == false) { //we are in the outer loop, global refinement is occuring

			t0 = t.back();

			//begin error bound calculation
			vector<bool> dumRef (dim, true);
			f(t0, y, dim, rowCount, k1, dumRef); //compute both components of Heun's
			for(int i = 0; i < dim; i++) k1[i] = k1[i] * h;
			vector<double> k2temp (dim, 0.0);
			for(int i = 0; i < dim; i++) k2temp[i] = y[i + (dim * rowCount)] + k1[i];
			f(t0 + h, k2temp, dim, 0, k2, dumRef);
			for(int i = 0; i < dim; i++) k2[i] = h * k2[i];

			
			for(int i = 0; i < dim; i++){
				double err = 0.5 * (k2[i] - k1[i]);
				err = fabs(err); //compute error comparing forward and Heun

				double tolC = relTol * fabs(y[i + (dim * rowCount)]) + absTol;
				double nErr = err / tolC;
				nErrVec[i] = nErr; //stores that node's normalised error
			}

			//check for the maximum error
			double nErrMax = 0.0;
			double nref_nErrMax = 0.0;
			for(int i = 0; i < dim; i++){
				if(nErrMax <= nErrVec[i]) nErrMax = nErrVec[i];
			}

			vector<bool> ref (dim, false); //create a default refinement vector where no node needs refining
			
			//now we check for refinement
			if(nErrMax > 1.0){ //if the worst error is outside tolerance then some refinement is necessary
				for(int i = 0; i < dim; i++){
					if(nErrVec[i] > (deltaRK*nErrMax)) ref[i] = true; //identifies nodes which need refinement
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
			else if(nref_nErrMax > 1.0){ //the not to be refined set has an error which is too large and the global step needs repeating with rescaled time-step
				fac = nu * pow(nref_nErrMax, -0.5);
				fac = max(hLo, fac);
				h = fac * h;
			}
			else{
				if(check > 0){ //this checks if there are components that need refinement		
					double Tin = t0 + h; //sets a new truncated integrating endpoint for refinement purposes
					fac = nu * pow(nErrMax, -0.5);
					fac = max(hLo, fac);
					double hin = fac * h; //rescaled h for refinement 

					//interpolation variables
					vector<double> y1nref (dim, 0.0);
					//the following lines equates to y1nref = y + 0.5 * (k1 + k2) from above
					for(int i = 0; i < dim; i++){
						if(ref[i] == 0){ //i.e not to be refined 
							y1nref[i] =  y[i + (dim * rowCount)] + 0.5 * (k1[i] + k2[i]);
						}
					}

					//assign information to be used by the inner loop

					inter.y1nref = y1nref;
					inter.t0 = t;
					inter.y0 = y;
					inter.h = h;
					inter.ref = ref;
					inter.OLDrowCount = rowCount;
					
					//we need internal refinements so set the check as true
					inter.internalCheck = true;
					
					//recursively call the function for refinement
					//the check guarantees we enter into the inner refinement stage
					solRK12MR(t, y, hin, Tin, f, relTol, absTol, dim, inter, rowCount);

					inter.internalCheck = false;

				}
				//else if the refinement vector is empty then we are good to proceed with appending the solutions
				else{
					//initialise an "empty" row onto the y vector (this might not even be necessary)
					for (int extend = 0; extend < dim; extend++) y.push_back(0.0);
				
					//assign the valid solutions to the y-vector
					for(int i = 0; i < dim; i++){ 
						y[i + (dim * (rowCount + 1))] = y[i + (dim * rowCount)] + 0.5 * (k1[i] + k2[i]);
					}

					//append t
					t.push_back(t0+h);
					t0 = t0 + h;
					fac = nuG * pow(nErrMax, -0.5);
					fac = min(hHiG, fac);
					h = h * fac;
					if(t0 + h > T) h = T - t0;									
					
					//iterate the rowCount variable to allow for smooth assigning onto the y-vector later
					rowCount++;

				}			
			}
		}
		//we are in the inner loop and need refinement
		else {
			t0 = t.back();
	
			//for nodes that need refinement compute both components of Heun's
			vector<double> k1 (dim, 0.0);
			f(t0, y, dim, rowCount, k1, inter.ref); //compute both components of Heun's
			for(int i = 0; i < dim; i++) k1[i] = k1[i] * h;

			vector<double> k2temp (dim, 0.0);
			for(int i = 0; i < dim; i++) k2temp[i] = (inter.ref[i] == true) ?  y[i + (dim * rowCount)] + k1[i] : inter.y0[i + (dim *inter.OLDrowCount)] + ((t0+h) - inter.t0.back())*(inter.y1nref[i] - inter.y0[i + (dim * inter.OLDrowCount)])/inter.h;
			f(t0 + h, k2temp, dim, 0, k2, inter.ref);
			for(int i = 0; i < dim; i++) k2[i] = h * k2[i];
		
			//error computation
			for(int i = 0; i < dim; i++){
				double err = 0.5 * (k2[i] - k1[i]);
				err = fabs(err); //compute error comparing forward and Heun's
				
				double tolC = relTol * fabs(y[i + (dim * rowCount)]) + absTol;
				double nErr = err / tolC;
				nErrVec[i] = nErr; //stores that node's normalised error
			}

			//check for the maximum error
			double nErrMax = 0.0;
			double nref_nErrMax = 0.0;
			for(int i = 0; i < dim; i++){
				if(nErrMax <= nErrVec[i]) nErrMax = nErrVec[i];
			}

			vector<bool> ref (dim, false); //create a default refinement vector where no node needs refining
			
			//now we check for refinement
			if(nErrMax > 1.0){ //if the worst error is outside tolerance then some refinement is necessary
				for(int i = 0; i < dim; i++){
					if(nErrVec[i] > (deltaRK*nErrMax)) ref[i] = true; //identifies nodes which need refinement
				}
			}
			
			//we now need to identify the max error of the NOT to be refined set
			for(int i = 0; i < dim; i++){
				if(ref[i] == false && inter.ref[i] == true){ //we only care about components that were marked for refinement in the outer loop
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
				//need to reset ref back to start of refinement process
				ref = inter.ref;
			}
			else if(nref_nErrMax > 1.0){ //nref error too large and the refinement process deemed to have failed
				fac = nu * pow(nref_nErrMax, -0.5);
				fac = max(hLo, fac);
				h = fac * h;
				//need to reset ref back to the start of the refinement process
				ref = inter.ref;
			}

			else{
				if(check > 0){ //this checks if there are components that need refinement

					double Tin = t0 + h;
					fac = nu * pow(nErrMax, -0.5);
					fac = max(hLo, fac);
					double hin = fac * h;

					//interpolation variables
					vector<double> y1nref (dim, 0.0);
					for(int i = 0; i < dim; i++){
						if(ref[i] == false && inter.ref[i] == true){ 
							y1nref[i] = y[i + (dim * rowCount)] + 0.5 * (k1[i] + k2[i]);

						}
						else if(inter.ref[i] == false){
							y1nref[i] += ((Tin) - inter.t0.back())*(inter.y1nref[i] - inter.y0[i + (dim * inter.OLDrowCount)])/inter.h;
						}
					}

					internalInfo inter1;

					inter1.y1nref = y1nref;
					inter1.t0 = t;
					inter1.y0 = y;
					inter1.h = h;
					inter1.ref = ref;
					inter1.OLDrowCount = rowCount;

					inter1.internalCheck = true;
					
					//we recursively call the function
					solRK12MR(t, y, hin, Tin, f, relTol, absTol, dim, inter1, rowCount);
					fac = nu * pow(nref_nErrMax , -0.5);
					fac = min(hHi, fac);
					h = fac * h;
					
					t0 = t.back();

					if(t0 + h > T) h = T - t0; 
				}
				else{
					//initialise an "empty" row onto the y vector
					ref = inter.ref;
					for (int extend = 0; extend < dim; extend++) y.push_back(0.0);
				
					for(int i = 0; i < dim; i++){
						if(ref[i] == 1){ 
							y[i + (dim * (rowCount + 1))] = y[i + (dim * rowCount)] + 0.5 * (k1[i] + k2[i]);
						}
						else{
							y[i + (dim * (rowCount + 1))] = inter.y0[i+(dim*inter.OLDrowCount)] + ((t0 + h) - inter.t0.back())*(inter.y1nref[i] - inter.y0[i + (dim * inter.OLDrowCount)])/inter.h;
						}
					}
					//append t
					t0 = t0 + h;
					t.push_back(t0);

					fac = nu * pow(nErrMax, -0.5);
					fac = min(hHi, fac);
					h = fac * h;

					if(t0 + h > T) h = T - t0;									
					
					rowCount++;
				}
			} 
		}
	}
}
