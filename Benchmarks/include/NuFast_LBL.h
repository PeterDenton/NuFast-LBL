#ifndef NuFast_LBL_H
#define NuFast_LBL_H

#include <vector>
#include <array>

using std::vector, std::array;

namespace NuFast {

void Probability_Matter_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, vector<double> E, double rhoYe, int N_Newton, vector<array<array<double, 3>, 3>> *probs_returned);
void Probability_Vacuum_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, vector<double> E, double empty1, int empty2, vector<array<array<double, 3>, 3>> *probs_returned);
}

#endif
