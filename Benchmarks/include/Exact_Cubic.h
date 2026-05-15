#ifndef Exact_Cubic_H
#define Exact_Cubic_H

#include <vector>
#include <array>

using std::vector, std::array;

void Probability_Matter_LBL_Exact_Cubic(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, const vector<double> &E, double rhoYe, int empty, vector<array<array<double, 3>, 3>> &probs_returned);

#endif
