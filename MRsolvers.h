#pragma once
#include<vector>;


//this is a very rough outline 
//will be substantially changed soon

struct internalInfo {
	bool internalCheck = false;
	float y1nref;
	float t0;
	float y0;
	float h;
	float ref;
	float tref;
};

void solRK12MR(std::vector<float>& t, std::vector<float>& y, float h, float T, float (*f)(float, float), internalInfo inter);