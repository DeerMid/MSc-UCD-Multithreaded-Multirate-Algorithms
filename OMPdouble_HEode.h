#include<vector>




void HEode(double t0, std::vector<double>& y, int dim,int row, std::vector<double>& ans, std::vector<bool>& ref);

void HEodeNOREF(double t0, std::vector<double>& y, int dim, int row, std::vector<double>& ans);
