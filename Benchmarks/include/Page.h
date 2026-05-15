#ifndef Page_H
#define Page_H

#include <vector>
#include <array>

using std::vector, std::array;

void Probability_Matter_LBL_Page(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, const vector<double> &E, double rhoYe, int empty, vector<array<array<double, 3>, 3>> &probs_returned);

#endif
