#include <cmath>
#include <vector>
#include <array>

#include "Page.h"
#include "Benchmark.h"

using std::vector, std::array;

#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))

void Probability_Matter_LBL_Page(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, const vector<double> &E, double rhoYe, int empty, vector<array<array<double, 3>, 3>> &probs_returned)
{

	// compute trig functions of mixing angles
	double c12sq, s12, c12, c13sq, s13, c13, c23sq, s23, c23, cd, sd;
	c12sq = 1 - s12sq;
	s12 = sqrt(s12sq);
	c12 = sqrt(c12sq);
	c13sq = 1 - s13sq;
	s13 = sqrt(s13sq);
	c13 = sqrt(c13sq);
	c23sq = 1 - s23sq;
	s23 = sqrt(s23sq);
	c23 = sqrt(c23sq);
	cd = cos(delta);
	sd = sin(delta);

	// Compute vacuum values first, independent of a, L or E
	double a0_vac, a1_vac, H_ee_tilde, Y_ee_tilde, R_H_em, I_H_em, R_Y_em, I_Y_em, PHASE;
	// Compute vacuum constants common to all flavours
	a0_vac = (Dmsq21*Dmsq21*Dmsq21 + Dmsq31*Dmsq31*Dmsq31) / 27.0 - (Dmsq21*Dmsq21 * Dmsq31 + Dmsq21 * Dmsq31*Dmsq31) / 18.0;
	a1_vac = (Dmsq21*Dmsq21 + Dmsq31*Dmsq31 - Dmsq21 * Dmsq31) / 9.0;
	H_ee_tilde = Dmsq21 * (s12*s12 * c13*c13 - 1.0 / 3.0) + Dmsq31 * (s13*s13 - 1.0 / 3.0);
	Y_ee_tilde = (Dmsq21*Dmsq21 * (s12*s12 * c13*c13 - 1.0 / 3.0) + Dmsq31*Dmsq31 * (s13*s13 - 1.0 / 3.0) + 2.0 * Dmsq21 * Dmsq31 * (c12*c12 * c13*c13 - 1.0 / 3.0)) / 3.0;

	// Compute extra mu->e constants
	R_H_em = Dmsq21 * s12 * c13 * (c12 * c23 - s12 * s23 * s13 * cd) + Dmsq31 * s13 * s23 * c13 * cd;
	I_H_em = Dmsq21 * s12*s12 * s13 * s23 * c13 * sd - Dmsq31 * s13 * s23 * c13 * sd;

	R_Y_em = (s12 * c13 * (c12 * c23 - s12 * s13 * s23 * cd) * Dmsq21*Dmsq21
			+ s13 * s23 * c13 * cd * Dmsq31*Dmsq31
			- 2.0 * c12 * c13 * (s12 * c23 + s13 * s23 * c12 * cd) * Dmsq21 * Dmsq31) / 3.0;

	I_Y_em = (s12*s12 * s13 * s23 * c13 * sd * Dmsq21*Dmsq21
			- s13 * s23 * c13 * sd * Dmsq31*Dmsq31
			+ 2.0 * s13 * s23 * c12*c12 * c13 * sd * Dmsq21 * Dmsq31) / 3.0;

	// Extra
	PHASE = 2.0 * M_PI / 3.0;

	// Compute all the matter-corrected quantities and then the transition probability
	double a, a0, a1, sqrt_a1, eigen[3], R_X[3], I_X[3], arcCos, L4E, Theta_10, Theta_20, Theta_21, denom;

	for (unsigned int i = 0; i < E.size(); i++)
	{
		a = rhoYe * E[i] * YerhoE2a;

		// Make  matter corrections (a0 and a1 have nothing to do with a)
		a /= 3.0;
		a0 = a0_vac + 1.5 * (Y_ee_tilde * a + H_ee_tilde * sq(a)) + cube(a); // a0
		a1 = a1_vac + H_ee_tilde * a + sq(a); // a1

		// Get eigenvalues of H, and constants X
		sqrt_a1 = sqrt(a1);
		arcCos = acos(a0 / (sqrt_a1 * a1)) / 3.0;
		sqrt_a1 *= 2.0;

		for (unsigned int j = 0; j < 3; ++j)
		{
			eigen[j] = sqrt_a1 * cos(arcCos - PHASE * j);
			denom = sq(eigen[j]) - a1;
			// pull out a 3 here, becomes a 9 later
			R_X[j] = ((eigen[j] + a) * R_H_em + R_Y_em) / denom;
			I_X[j] = ((eigen[j] + a) * I_H_em + I_Y_em) / denom;
		}

		// Compute Theta constants (not mixing angles)
		L4E = eVsqkm_to_GeV_over4 * L / E[i];
		Theta_10 = (eigen[1] - eigen[0]) * L4E;
		Theta_20 = (eigen[2] - eigen[0]) * L4E;
		Theta_21 = (eigen[2] - eigen[1]) * L4E;

		// Compute probabilities
		double Pme_CPC, Pme_CPV;
		// mu->e
		Pme_CPC = 2.0 * (- 2.0 * ((R_X[1]*R_X[0] + I_X[1]*I_X[0]) * sq(sin(Theta_10))
								+ (R_X[2]*R_X[0] + I_X[2]*I_X[0]) * sq(sin(Theta_20))
								+ (R_X[2]*R_X[1] + I_X[2]*I_X[1]) * sq(sin(Theta_21)))) / 9.0;
		Pme_CPV = 2.0 * (((I_X[1]*R_X[0] - R_X[1]*I_X[0]) * sin(2.0 * Theta_10)
								+ (I_X[2]*R_X[0] - R_X[2]*I_X[0]) * sin(2.0 * Theta_20)
								+ (I_X[2]*R_X[1] - R_X[2]*I_X[1]) * sin(2.0 * Theta_21))) / 9.0;

		// Compute the extra e->e constants
		// e->e
		double H_ee, Y_ee, X_ee[3];
		H_ee = H_ee_tilde + 2 * a; // note that "a" already has a 1/3 on it
		Y_ee = Y_ee_tilde + 2 * H_ee_tilde * a + 2 * sq(a);
		for (int j = 0; j < 3; j++)
			X_ee[j] = 1 + (eigen[j] * H_ee + Y_ee) / (sq(eigen[j]) - a1); // pull out a 3 here, becomes a 9 later
		probs_returned[i][0][0] = 1. - 4 * (X_ee[1] * X_ee[0] * sq(sin(Theta_10))
										  + X_ee[2] * X_ee[0] * sq(sin(Theta_20))
										  + X_ee[2] * X_ee[1] * sq(sin(Theta_21))) / 9.;

		probs_returned[i][1][0] = Pme_CPC + Pme_CPV;

		// mu->mu
		double H_mm_tilde, H_mm, Y_mm_tilde, Y_mm, X_mm[3];
		H_mm_tilde = Dmsq21 * (c12sq * c23sq + s12sq * s13sq * s23sq - 2 * s12 * s13 * s23 * c12 * c23 * cd - 1. / 3) + Dmsq31 * (s23sq * c13sq - 1. / 3);
		H_mm = H_mm_tilde - a;
		Y_mm_tilde = (sq(Dmsq21) * (c12sq * c23sq + s12sq * s13sq * s23sq - 2 * s12 * s13 * s23 * c12 * c23 * cd - 1. / 3) + sq(Dmsq31) * (s23sq * c13sq - 1. / 3) + 2 * Dmsq21 * Dmsq31 * (s12sq * c23sq + s13sq * s23sq * c12sq + 2 * s12 * s13 * s23 * c12 * c23 * cd - 1. / 3)) / 3;
		Y_mm = Y_mm_tilde - 2 * (H_ee_tilde + H_mm_tilde) * a - sq(a);
		for (int j = 0; j < 3; j++)
			X_mm[j] = 1 + (eigen[j] * H_mm + Y_mm) / (sq(eigen[j]) - a1); // pull out a 3 here, becomes a 9 later
		probs_returned[i][1][1] = 1. - 4 * (X_mm[1] * X_mm[0] * sq(sin(Theta_10))
										  + X_mm[2] * X_mm[0] * sq(sin(Theta_20))
										  + X_mm[2] * X_mm[1] * sq(sin(Theta_21))) / 9.;


		probs_returned[i][2][0] = 1. - probs_returned[i][1][0] - probs_returned[i][0][0]; // tau->e
		probs_returned[i][0][1] = Pme_CPC - Pme_CPV; // e->mu
		probs_returned[i][0][2] = 1. - probs_returned[i][0][0] - probs_returned[i][0][1]; // e->tau
		probs_returned[i][1][2] = 1. - probs_returned[i][1][0] - probs_returned[i][1][1]; // mu->tau
		probs_returned[i][2][1] = 1. - probs_returned[i][0][1] - probs_returned[i][1][1]; // tau->mu
		probs_returned[i][2][2] = 1. - probs_returned[i][2][0] - probs_returned[i][2][1]; // tau->tau
	} // i, energies
}

