#ifndef DMP_H
#define DMP_H

#include <vector>
#include <array>

using std::vector, std::array;

// From http://arxiv.org/abs/1604.08167
// Zeorth order only

void Probability_Matter_LBL_DMP(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, vector<double> E, double rhoYe, int empty, vector<array<array<double, 3>, 3>> *probs_returned);

#endif
