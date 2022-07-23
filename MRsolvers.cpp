#include "MRsolvers.h"
#include<vector>
#include<iostream>

using namespace std;

void solRK12MR(vector<float>& t, vector<float>& y, float h, float T, float (*f)(float, float), int dim, internalInfo inter) {
	
	float t0 = t.back();

	while (t0 < T) {
		if (inter.internalCheck == false) { //we are in the outer loop


		}
		else {

		}
	}
}
