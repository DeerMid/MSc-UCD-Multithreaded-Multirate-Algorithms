#pragma once
#include<vector>


//this is a very rough outline 
//will be substantially changed soon

struct internalInfo {
	bool internalCheck = false;
	std::vector<double> y1nref;
	std::vector<double> t0;
	std::vector<double> y0;
	double h;
	std::vector<bool> ref;
	int OLDrowCount = 0;
};

void solRK12MR(std::vector<double>& t, std::vector<double>& y, double h, double T, void (*f)(double, std::vector<double>&, int, int, std::vector<double>&, std::vector<bool>&), double relTol, double absTol, int dim, internalInfo &inter, int &rowCount);

