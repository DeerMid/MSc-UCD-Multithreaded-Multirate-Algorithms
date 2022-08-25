#include<stdio.h>
#include<vector>
#include<math.h>
#include<iostream>

using namespace std;

void HEode(double t0, vector<double>& y, int dim,int row, vector<double>& ans, vector<bool>& ref){

	const double k = 1.0e-2; 

	
	vector<double> yCopy (dim, 0.0);
	for(int i = 0; i < dim; i++) yCopy[i] = y[i + (dim * row)];

	//padding the boundary values
	yCopy.insert(yCopy.begin(), 0.0);
	yCopy.push_back(0.0);

	double dx = 1.0/dim;

	for(int i = 0; i < dim; i++){
		ans[i] = (ref[i] == true) ?  k * (yCopy[i+2] - 2*yCopy[i+1] + yCopy[i])/(pow(dx, 2)) : 0.0;
	}
	//for(int i = 0; i < dim; i++) cout << " HE ans[" << i << "]: " << ans[i];
	//cout << endl;

}

void HEodeNOREF(double t0, vector<double>& y, int dim, int row, vector<double>& ans){

	const double k = 1.0e-2; 

	vector<double> yCopy (dim, 0.0);
	for(int i = 0; i < dim; i++) yCopy[i] = y[i + (dim * row)];

	//padding the boundary values
	yCopy.insert(yCopy.begin(), 0.0);
	yCopy.push_back(0.0);

	double dx = 1.0/dim;

	for(int i = 0; i < dim; i++){
		ans[i] = yCopy[i+2] - 2*yCopy[i+1] + yCopy[i]/(pow(dx, 2));
	}

}
