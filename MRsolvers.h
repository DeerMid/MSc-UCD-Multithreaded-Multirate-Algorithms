#pragma once
#include<vector>


//this is a very rough outline 
//will be substantially changed soon

struct internalInfo {
	bool internalCheck = false;
	std::vector<float> y1nref;
	std::vector<float> t0;
	std::vector<float> y0;
	float h;
	std::vector<bool> ref;
	int rowCount = 0;
};

void solRK12MR(std::vector<float>& t, std::vector<float>& y, float h, float T, float (*f)(float, float), float relTol, float absTol, int dim, internalInfo inter);
