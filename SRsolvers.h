#pragma once
#include<vector>

void solForwardEulerSR(std::vector<float>& t, std::vector<float>& y, float h, float T, float (*f)(float, float), int dim);
void solHeunSR(std::vector<float>& t, std::vector<float>& y, float h, float T, float (*f)(float, float), int dim);
void solRK12SR(std::vector<float>& t, std::vector<float>& y, float h, float T, float (*f)(float, float), float relTol, float absTol, int dim);
