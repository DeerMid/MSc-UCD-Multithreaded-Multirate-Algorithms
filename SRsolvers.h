#pragma once
#include<vector>

void solForwardEulerSR(std::vector<float>& t, std::vector<float>& y, float h, float T, float (*f)(float, float), int dim);
void solHeunSR(std::vector<float>& t, std::vector<float>& y, float h, float T, float (*f)(float, float), int dim);
void solRK12SR(std::vector<double>& t, std::vector<double>& y, double h, double T, void (*f)(double, std::vector<double>&, int, int, std::vector<double>&, std::vector<bool>&), double relTol, double absTol, int dim);
