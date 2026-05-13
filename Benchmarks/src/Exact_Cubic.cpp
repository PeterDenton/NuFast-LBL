#include <cmath>
#include <vector>
#include <array>

#include "Exact_Cubic.h"
#include "Benchmark.h"

using std::vector, std::array;

void Probability_Matter_LBL_Exact_Cubic(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, vector<double> E, double rhoYe, int empty, vector<array<array<double, 3>, 3>> *probs_returned)
{
    double Amatter, c13sq, sind, cosd;

    double A, A_tmp, B, C, S;
    double See, Smm, Tee, Tmm, Tmm_tmp, rootAsqm3B, tmp;
    double lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32;
    double Ue1sq, Ue2sq, Ue2sq_tmp, Ue3sq, Ue3sq_tmp, Um1sq, Um2sq, Um2sq_tmp, Um3sq, Um3sq_tmp, Ut1sq, Ut2sq, Ut2sq_tmp, Ut3sq;

    double sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin, PiDlambdaInv, Xp2, Xp3;
    double D21, D32;
    double sinD21, sinD31, sinD32;
    double Lover4E, Jrr, Jmatter, Jmatter_tmp;
    double Pee, Pme_CPC, Pme_CPV, Pmm;

    // --------------------------------------------------------------------- //
    // First calculate useful simple functions of the oscillation parameters //
    // --------------------------------------------------------------------- //
    c13sq = 1 - s13sq;

    // Ueisq's
    Ue2sq_tmp = c13sq * s12sq;
    Ue3sq_tmp = s13sq;

    // Umisq's, Utisq's and Jvac
    Um3sq_tmp = c13sq * s23sq;
    Ut2sq_tmp = s13sq * s12sq * s23sq;
    Um2sq_tmp = (1 - s12sq) * (1 - s23sq);
      
    Jrr = sqrt(Um2sq_tmp * Ut2sq_tmp);
    sind = sin(delta);
    cosd = cos(delta);

    Um2sq_tmp = Um2sq_tmp + Ut2sq_tmp - 2 * Jrr * cosd;
    Jmatter_tmp = 8 * Jrr * c13sq * sind;

	A_tmp = Dmsq21 + Dmsq31; // temporary variable
	See = A_tmp - Dmsq21 * Ue2sq_tmp - Dmsq31 * Ue3sq_tmp;
	Tmm_tmp = Dmsq21 * Dmsq31; // using Tmm as a temporary variable      
	Tee = Tmm_tmp * (1 - Ue3sq_tmp - Ue2sq_tmp);

	for (unsigned int i = 0; i < E.size(); i++)
	{
		Amatter = rhoYe * E[i] * YerhoE2a;

		C = Amatter * Tee;
		A = A_tmp + Amatter;

		// -------------------------------- //
		// Get exact lambda3 from the cubic //
		// -------------------------------- //
		B = Tmm_tmp + Amatter * See;
		rootAsqm3B = sqrt(A * A - 3 * B);
		S = acos((A * A * A - 4.5 * A* B + 13.5 * C) / (rootAsqm3B * rootAsqm3B * rootAsqm3B));
		if (Dmsq31 < 0) S = S + 2 * M_PI; // adjust lambda3 for the IO
		lambda3 = (A + 2 * rootAsqm3B * cos(S * oneThird)) * oneThird;

		// ------------------- //
		// Get  Delta lambda's //
		// ------------------- //
		tmp = A - lambda3;
		Dlambda21 = sqrt(tmp * tmp - 4 * C / lambda3);
		lambda2 = 0.5 * (A - lambda3 + Dlambda21);
		Dlambda32 = lambda3 - lambda2;
		Dlambda31 = Dlambda32 + Dlambda21;

		// ----------------------- //
		// Use Rosetta for Ueisq's //
		// ----------------------- //
		// denominators      
		PiDlambdaInv = 1 / (Dlambda31 * Dlambda32 * Dlambda21);
		Xp3 = PiDlambdaInv * Dlambda21;
		Xp2 = -PiDlambdaInv * Dlambda31;

		// numerators
		Ue3sq = (lambda3 * (lambda3 - See) + Tee) * Xp3;
		Ue2sq = (lambda2 * (lambda2 - See) + Tee) * Xp2;

		Smm = A - Dmsq21 * Um2sq_tmp - Dmsq31 * Um3sq_tmp;
		Tmm = Tmm_tmp * (1 - Um3sq_tmp - Um2sq_tmp) + Amatter * (See + Smm - A);

		Um3sq = (lambda3 * (lambda3 - Smm) + Tmm) * Xp3;
		Um2sq = (lambda2 * (lambda2 - Smm) + Tmm) * Xp2;

		// ------------- //
		// Use NHS for J //
		// ------------- //
		Jmatter = Jmatter_tmp * Dmsq21 * Dmsq31 * (Dmsq31 - Dmsq21) * PiDlambdaInv;

		// ----------------------- //
		// Get all elements of Usq //
		// ----------------------- //
		Ue1sq = 1 - Ue3sq - Ue2sq;
		Um1sq = 1 - Um3sq - Um2sq;

		Ut3sq = 1 - Um3sq - Ue3sq;
		Ut2sq = 1 - Um2sq - Ue2sq;
		Ut1sq = 1 - Um1sq - Ue1sq;

		// ----------------------- //
		// Get the kinematic terms //
		// ----------------------- //
		Lover4E = eVsqkm_to_GeV_over4 * L / E[i];

		D21 = Dlambda21 * Lover4E;
		D32 = Dlambda32 * Lover4E;
		  
		sinD21 = sin(D21);
		sinD31 = sin(D32 + D21);
		sinD32 = sin(D32);

		triple_sin = sinD21 * sinD31 * sinD32;

		sinsqD21_2 = 2 * sinD21 * sinD21;
		sinsqD31_2 = 2 * sinD31 * sinD31;
		sinsqD32_2 = 2 * sinD32 * sinD32;

		// ------------------------------------------------------------------- //
		// Calculate the three necessary probabilities, separating CPC and CPV //
		// ------------------------------------------------------------------- //
		Pme_CPC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
				+ (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
				+ (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
		Pme_CPV = -Jmatter * triple_sin;

		Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2
					 + Um3sq * Um1sq * sinsqD31_2
					 + Um3sq * Um2sq * sinsqD32_2);

		Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2
					 + Ue3sq * Ue1sq * sinsqD31_2
					 + Ue3sq * Ue2sq * sinsqD32_2);

		// ---------------------------- //
		// Assign all the probabilities //
		// ---------------------------- //
		(*probs_returned)[i][0][0] = Pee;															// Pee
		(*probs_returned)[i][0][1] = Pme_CPC - Pme_CPV;												// Pem
		(*probs_returned)[i][0][2] = 1 - Pee - (*probs_returned)[i][0][1];  						// Pet

		(*probs_returned)[i][1][0] = Pme_CPC + Pme_CPV;												// Pme
		(*probs_returned)[i][1][1] = Pmm;															// Pmm
		(*probs_returned)[i][1][2] = 1 - (*probs_returned)[i][1][0] - Pmm;							// Pmt

		(*probs_returned)[i][2][0] = 1 - Pee - (*probs_returned)[i][1][0];							// Pte
		(*probs_returned)[i][2][1] = 1 - (*probs_returned)[i][0][1] - Pmm;							// Ptm
		(*probs_returned)[i][2][2] = 1 - (*probs_returned)[i][0][2] - (*probs_returned)[i][1][2];	// Ptt
	} // i, energies
}
