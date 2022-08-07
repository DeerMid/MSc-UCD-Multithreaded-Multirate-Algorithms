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
	float nErrMax = 1.0; //preinitialise
	float nref_nErrMax = 0.0;

	//create vectors to contain the components of Heun's method
	vector<float> k1;
	vector<float> k2;

	vector<float> nErrVec; //vector to store a time step's normalised error	

	float t0 = t.back();

	while (t0 < T) {
		if (inter.internalCheck == false) { //we are in the outer loop
		
			float t0 = t.back();
	
			//the y-push back will most likely need to be delayed until actual assignment later, otherwise y expands a row every function call!
			//although if push back absent in internal refinement then this should be fine
			//will confirm later!
			//for (int extend = 0; extend < dim; extend++) y.push_back(0.0);
			for (int i = 0; i < dim; i++) {
				k1[i] = h * f(t0, y[i + (dim * inter.rowCount)]); //compute both components of Heun's
				k2[i] = h * f(t0 + h, y[i + (dim * inter.rowCount)] + k1[i]);

				float err = 0.5 * (k1[i] - k2[i]);
				err = fabs(err); //compute error comparing forward and Heun's

				float tolC = relTol * fabs(y[i + (dim * inter.rowCount)]) + absTol;
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
					vector<float> yInter (dim, 0.0);
					//the following line equates to yInter = y + 0.5 * (k1 + k2) from above
					for(int i = 0; i < dim; i++) yInter[i] = y[i + (dim * inter.rowCount)] + 0.5 * (k1[i] + k2[i]);
					vector<float> y1nref;
					for(int i = 0; i < dim; i++){
						if(ref[i] == 0){ 
							y1nref[i] = yInter[i] * 1.0;
						}
					}

					//tref may be needed here

					inter.y1nref = y1nref;
					inter.t0 = t;
					inter.y0 = y;
					inter.h = h;
					inter.ref = ref;

					inter.internalCheck = true;
					
					//recursively call the function for refinement
					solRK12MR(t, y, hin, Tin, f, relTol, absTol, dim, inter);  
				}
				else{
					//initialise an "empty" row onto the y vector
					for (int extend = 0; extend < dim; extend++) y.push_back(0.0);
				
					for(int i = 0; i < dim; i++){ 
						y[i + (dim * (inter.rowCount + 1))] = y[i + (dim * inter.rowCount)] + 0.5 * (k1[i] + k2[i]);
					}

					//append t
					t0 = t0 + h;
					t.push_back(t0);									
					
					inter.rowCount++;
				}			
			}

			

			//debugging tool, force closes the while loop and prints a message 
			cout << "Made it this far!" << endl;
			t0 = T + 100;
		}
		else {
		
			float t0 = t.back();
	
			//unpack the structure
			//vector<float> y1r = inter.y1nref;
			//vector<float> t0r = inter.t0;
			//vector<float> y0r = inter.y0;
			//float hr = inter.h;
			//vector<bool> ref
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
					vector<float> yInter (dim, 0.0);
					//the following line equates to yInter = y + 0.5 * (k1 + k2) from above
					for(int i = 0; i < dim; i++) yInter[i] = y[i + (dim * inter.rowCount)] + 0.5 * (k1[i] + k2[i]);
					vector<float> y1nref;
					for(int i = 0; i < dim; i++){
						if(ref[i] == 0){ 
							y1nref[i] = yInter[i] * 1.0;
						}
						else{
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
					solRK12MR(t, y, hin, Tin, f, relTol, absTol, dim, inter);

					fac = nu * pow(nref_nErrMax , -0.5);
					fac = min(hHi, fac);
					h = fac * h;

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
			cout << "This is the internal refinement stage, come back later!" << endl;
		}
	}
}
