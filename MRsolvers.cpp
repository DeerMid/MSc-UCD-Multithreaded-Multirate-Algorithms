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
	const float nuG = 1.0;
	const float hHiG = 1.0;
	float fac = 1.0; 

	float nErrMax = 1.0; //preinitialise error maximums
	float nref_nErrMax = 0.0;

	//create vectors to contain the components of Heun's method
	vector<float> k1 (dim, 0.0);
	vector<float> k2 (dim, 0.0);

	vector<float> nErrVec (dim, 0.0); //vector to store a time step's normalised error	

	float t0 = t.back();

	while (t0 < T) {
		if (inter.internalCheck == false) { //we are in the outer loop, global refinement is occuring
			cout << "Entered outer loop" << endl;	
			float t0 = t.back();
			cout << "t0 is " << t0 << endl;
			//begin error bound calculation
			for (int i = 0; i < dim; i++) {
				//cout << i << " k1 and k2 calulation " << endl;
				k1[i] = h * f(t0, y[i + (dim * inter.rowCount)]); //compute both components of Heun's
				k2[i] = h * f(t0 + h, y[i + (dim * inter.rowCount)] + k1[i]);

				float err = 0.5 * (k1[i] - k2[i]);
				err = fabs(err); //compute error comparing forward and Heun's
				cout << "the error after fabs is " << err << endl;

				float tolC = relTol * fabs(y[i + (dim * inter.rowCount)]) + absTol;
				float nErr = err / tolC;
				nErrVec[i] = nErr; //stores that node's normalised error
				cout << nErrVec[i] << " is the error here" << endl;
			}
			//cout << "Completed calculations of k1 and k2 in outer loop" << endl;
			//check for the maximum error
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
			cout << "nref_nErrMax is " << nref_nErrMax << endl;
			
			//this step determines how "large" the refinement will be, this can probably be improved much further in terms of efficiency
			int check = 0;
			for(int i = 0; i < dim; i++) check += ref[i];
			cout << "the check for this loop is " << check << endl;

			if(check == dim){ //all components need to be refined
				fac = nu * pow(nErrMax, -0.5);
				fac = max(hLo, fac);
				h = fac * h;
				cout << "all components need refining" << endl;
			}
			else if(nref_nErrMax > 1.0){ //the not to be refined set has an error which is too large and the global step needs repeating with rescaled time-step
				fac = nu * pow(nref_nErrMax, -0.5);
				fac = max(hLo, fac);
				h = fac * h;
				cout << "global step failed" << endl;
			}
			else{
				if(check > 0){ //this checks if there are components that need refinement
					cout << "Refinement recursion calling" << endl;		
					float Tin = t0 + h; //sets a new truncated integrating endpoint for refinement purposes
					fac = nu * pow(nErrMax, -0.5);
					fac = max(hLo, fac);
					float hin = fac * h; //rescaled h for refinement 

					//interpolation variables
					vector<float> y1nref (dim, 0.0);
					//the following lines equates to y1nref = y + 0.5 * (k1 + k2) from above
					for(int i = 0; i < dim; i++){
						if(ref[i] == 0){ //i.e not to be refined 
							y1nref[i] =  y[i + (dim * inter.rowCount)] + 0.5 * (k1[i] + k2[i]);;
						}
					}

					//assign information to be used by the inner loop

					inter.y1nref = y1nref;
					inter.t0 = t;
					inter.y0 = y;
					inter.h = h;
					inter.ref = ref;
					
					//we need internal refinements so set the check as true
					inter.internalCheck = true;
					
					//recursively call the function for refinement
					//the check guarantees we enter into the inner refinement stage
					cout << "recursively calling function" << endl;
					solRK12MR(t, y, hin, Tin, f, relTol, absTol, dim, inter); 
					cout << "finished outer recursion safely" << endl;
					t0 = t.back();

					if(t0 + h > T) h = T - t0; 

				}
				//else if the refinement vector is empty then we are good to proceed with appending the solutions
				else{
					cout << "outer loop okay, appending and continuing" << endl;
					//initialise an "empty" row onto the y vector (this might not even be necessary)
					for (int extend = 0; extend < dim; extend++) y.push_back(0.0);
				
					//assign the valid solutions to the y-vector
					for(int i = 0; i < dim; i++){ 
						y[i + (dim * (inter.rowCount + 1))] = y[i + (dim * inter.rowCount)] + 0.5 * (k1[i] + k2[i]);
					}

					//append t
					t0 = t0 + h;
					t.push_back(t0);									
					
					//iterate the rowCount variable to allow for smooth assigning onto the y-vector later
					inter.rowCount++;

					cout << "t0 is, at the end, " << t0 << endl;
				}			
			}

		}
		//we are in the inner loop and need refinement
		else {
			cout << "in the inner refinement loop" << endl;
			float t0 = t.back();
	
			for(int i = 0; i < dim; i++){
				if(inter.ref[i] == true){
					//for nodes that need refinement compute both components of Heun's
					k1[i] = h * f(t0, y[i + (dim * inter.rowCount)]);
					k2[i] = h * f(t0 + h, y[i + (dim * inter.rowCount)] + k1[i]);
					
					//error computation
					float err = 0.5 * (k1[i] - k2[i]);
					err = fabs(err); //compute error comparing forward and Heun's
			
					float tolC = relTol * fabs(y[i + (dim * inter.rowCount)]) + absTol;
					float nErr = err / tolC;
					nErrVec[i] = nErr; //stores that node's normalised error
				}
				else{
					//interpolation for the nodes that don't need refinement and computing Heun's
					k1[i] = inter.y0[i] + ((t0+h) - inter.t0.back())*(inter.y1nref[i] - inter.y0[i])/inter.h;
					k2[i] = h * f(t0 + h, k1[i]);
					
					//error is irrelevant for these components
					nErrVec[i] = 0.0;  						
				}
			}
			//check for the maximum error
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

					float Tin = t0 + h;
					fac = nu * pow(nErrMax, -0.5);
					fac = max(hLo, fac);
					float hin = fac * h;

					//interpolation variables
					vector<float> yInter (dim, 0.0);
					//the following line equates to yInter = y + 0.5 * (k1 + k2) from above
					for(int i = 0; i < dim; i++) yInter[i] = y[i + (dim * inter.rowCount)] + 0.5 * (k1[i] + k2[i]);
					vector<float> y1nref (dim, 0.0);
					for(int i = 0; i < dim; i++){
						if(ref[i] == false && inter.ref[i] == true){ 
							y1nref[i] = y[i + (dim * inter.rowCount)] + 0.5 * (k1[i] + k2[i]);

						}
						else if(inter.ref[i] == true){
							y1nref[i] = ((Tin) - inter.t0.back())*(inter.y1nref[i] - inter.y0[i])/inter.h;
						}
					}
					

					//tref may be needed here

					inter.y1nref = y1nref;
					inter.t0 = t;
					inter.y0 = y;
					inter.h = h;
					inter.ref = ref;

					inter.internalCheck = true;
					
					//we recursively call the function
					cout << "calling recursively on inner loop" << endl;
					solRK12MR(t, y, hin, Tin, f, relTol, absTol, dim, inter);
					cout << "finished inner recursion" << endl;
					fac = nu * pow(nref_nErrMax , -0.5);
					fac = min(hHi, fac);
					h = fac * h;
					
					t0 = t.back();

					if(t0 + h > T) h = T - t0; 
				}
				else{
					//initialise an "empty" row onto the y vector
					for (int extend = 0; extend < dim; extend++) y.push_back(0.0);
				
					for(int i = 0; i < dim; i++){
						if(ref[i] == 0){ 
							y[i + (dim * (inter.rowCount + 1))] = y[i + (dim * inter.rowCount)] + 0.5 + (k1[i] + k2[i]);
						}
						else{
							y[i + (dim * (inter.rowCount + 1))] = ((t0 + h) - inter.t0.back())*(inter.y1nref[i] - inter.y0[i])/inter.h;
						}
					}

					//append t
					t0 = t0 + h;
					t.push_back(t0);

					fac = nuG * pow(nErrMax, -0.5);
					fac = min(hHiG, fac);
					h = fac * h;

					if(t0 + h > T) h = T - t0;									
					
					inter.rowCount++;
					inter.internalCheck = false; //set it so we continue in the outer refinement loop
				}
			}
		}
	}
}
