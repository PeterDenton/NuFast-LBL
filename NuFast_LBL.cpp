#include <cmath>
#include <stdio.h>
#include <vector>
#include <array>

using std::vector, std::array;

namespace NuFast {

// Some constants
constexpr double const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4;
constexpr double const YerhoE2a = 1.52588e-4;

// Probability_Matter_LBL calculates all nine oscillation probabilities across a
// vector of energies, including the matter effect, in an optimized, fast, and
// efficient way. The precision can be controlled with N_Newton. For many
// applications N_Newton=0 will be enough, but many years of DUNE or HK-LBL may
// require N_Newton=1 for a high precision nuerical scan. The code is standalone.
//
// Inputs:
//   mixing angles (usual parameterization)
//   phase (usual parameterization, radians)
//   Delta msq's (eV^2)
//   Make Dmsq31 positive/negative for the NO/IO
//   L (km)
//   E (GeV) positive for neutrinos, negative for antineutrinos. Vector
//   rhoYe density (g/cc) times electron fraction, typically around 3 * 0.5 for the Earth's crust
//   N_Newton: number of Newton's method iterations to do. Should be zero, one, two (can be higher); if it is negative, the code will use the exact expression
// Outputs:
//   probs_returned: a vector (over energies) of all nine oscillation probabilities: e.g. probs_returned[i][1][0] is P(nu_mu->nu_e) for energy E[i]
//     the vector should be allocated first
void Probability_Matter_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, const vector<double> &E, double rhoYe, int N_Newton, vector<array<array<double, 3>, 3>> &probs_returned)
{
	// ------------------------ //
	// The variables to be used //
	// ------------------------ //
	double c13sq, sind, cosd, Jrr, Jmatter, Jmatter_tmp, Dmsqee, Amatter;
	double Ue1sq, Ue2sq, Ue2sq_tmp, Ue3sq, Ue3sq_tmp, Um1sq, Um2sq, Um2sq_tmp, Um3sq, Um3sq_tmp,Ut1sq, Ut2sq, Ut2sq_tmp, Ut3sq;
	double A, A_tmp, B, C;
	double See, Tee, Smm, Tmm, Tmm_tmp;
	double rootAsqB, ss0;
	double xmat, lambda2, lambda3, Dlambda21, Dlambda31, Dlambda32;
	double Xp2, Xp3, PiDlambdaInv, tmp;
	double Lover4E, D21, D32;
	double sinD21, sinD31, sinD32;
	double sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin;
	double Pme_TC, Pme_TV, Pmm, Pee;

	// ------------------------------------------------------------------------------------------------- //
	// First calculate useful simple functions of the oscillation parameters that don't depend on energy //
	// ------------------------------------------------------------------------------------------------- //
	c13sq = 1 - s13sq;

	// Ueisq's
	Ue2sq_tmp = c13sq * s12sq;
	Ue3sq_tmp = s13sq;

	// Umisq's, Utisq's and Jvac
	Um3sq_tmp = c13sq * s23sq;
	// Um2sq_tmp and Ut2sq_tmp are used here as temporary variables, will be properly defined later	 
	Ut2sq_tmp = s13sq * s12sq * s23sq;
	Um2sq_tmp = (1 - s12sq) * (1 - s23sq);

	Jrr = sqrt(Um2sq_tmp * Ut2sq_tmp);
	sind = sin(delta);
	cosd = cos(delta);

	Um2sq_tmp += Ut2sq_tmp - 2 * Jrr * cosd;
	Jmatter_tmp = 8 * Jrr * c13sq * sind;
	Dmsqee = Dmsq31 - s12sq * Dmsq21;

	// calculate A, B, C, See, Tee, and part of Tmm
	A_tmp = Dmsq21 + Dmsq31; // temporary variable
	See = A_tmp - Dmsq21 * Ue2sq_tmp - Dmsq31 * Ue3sq_tmp;
	Tmm_tmp = Dmsq21 * Dmsq31; // using Tmm as a temporary variable	  
	Tee = Tmm_tmp * (1 - Ue3sq_tmp - Ue2sq_tmp);

	// --------------------------- //
	// Loop over the energy vector //
	// --------------------------- //
	for (unsigned int i = 0; i < E.size(); i++)
	{
		Amatter = rhoYe * E[i] * YerhoE2a;
		C = Amatter * Tee;
		A = A_tmp + Amatter;
		B = Tmm_tmp + Amatter * See; // B is not needed for N_Newton=0, but we calculate it everytime anyway

		if (N_Newton < 0)
		{
			// --------------------------------------------------------- //
			// Get lambda3 from exact expression, computationally costly //
			// Lambda3 for both mass orderings                           //
			// --------------------------------------------------------- //
			rootAsqB = sqrt(A * A - 3. * B);
			ss0 = acos((A * A * A - 4.5 * A * B + 13.5 * C) / (rootAsqB * rootAsqB * rootAsqB));
			if (Dmsq31 < 0.)
				ss0 = ss0 + 2. * M_PI; //  add 2Pi if Dmsq31 < 0 to get lambda3
			lambda3 = (A + 2. * rootAsqB * cos(ss0 / 3.)) / 3.;
		}
		else
		{
			// ---------------------------------- //
			// Get lambda3 from lambda+ of MP/DMP //
			// ---------------------------------- //
			xmat = Amatter / Dmsqee;
			tmp = 1 - xmat;
			lambda3 = Dmsq31 + 0.5 * Dmsqee * (xmat - 1 + sqrt(tmp * tmp + 4 * s13sq * xmat));

			// ---------------------------------------------------------------------------- //
			// Newton iterations to improve lambda3 arbitrarily, if needed, (B needed here) //
			// ---------------------------------------------------------------------------- //
			for (int j = 0; j < N_Newton; j++)
				lambda3 = (lambda3 * lambda3 * (lambda3 + lambda3 - A) + C) / (lambda3 * (2 * (lambda3 - A) + lambda3) + B); // this strange form prefers additions to multiplications
		}

		// ------------------- //
		// Get  Delta lambda's //
		// ------------------- //
		tmp = A - lambda3;
		Dlambda21 = sqrt(tmp * tmp - 4 * C / lambda3);
		lambda2 = 0.5 * (A - lambda3 + Dlambda21);
		Dlambda32 = lambda3 - lambda2;
		Dlambda31 = Dlambda32 + Dlambda21;

		// ----------------------- //
		// Use Rosetta for Veisq's //
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

		// --------------------------------- //
		// Get all elements of Usq in matter //
		// --------------------------------- //
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

		// ----------------------------------------------------------------- //
		// Calculate the three necessary probabilities, separating TC and TV //
		// ----------------------------------------------------------------- //
		Pme_TC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
			   + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
			   + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
		Pme_TV = -Jmatter * triple_sin;

		Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2
					 + Um3sq * Um1sq * sinsqD31_2
					 + Um3sq * Um2sq * sinsqD32_2);

		Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2
					 + Ue3sq * Ue1sq * sinsqD31_2
					 + Ue3sq * Ue2sq * sinsqD32_2);

		// ---------------------------- //
		// Assign all the probabilities //
		// ---------------------------- //
		probs_returned[i][0][0] = Pee;														// Pee
		probs_returned[i][0][1] = Pme_TC - Pme_TV;											// Pem
		probs_returned[i][0][2] = 1 - Pee - probs_returned[i][0][1];  						// Pet

		probs_returned[i][1][0] = Pme_TC + Pme_TV;											// Pme
		probs_returned[i][1][1] = Pmm;														// Pmm
		probs_returned[i][1][2] = 1 - probs_returned[i][1][0] - Pmm;						// Pmt

		probs_returned[i][2][0] = 1 - Pee - probs_returned[i][1][0];						// Pte
		probs_returned[i][2][1] = 1 - probs_returned[i][0][1] - Pmm;						// Ptm
		probs_returned[i][2][2] = 1 - probs_returned[i][0][2] - probs_returned[i][1][2];	// Ptt
	} // i, energies
}

void Probability_Vacuum_LBL(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, const vector<double> &E, vector<array<array<double, 3>, 3>> &probs_returned)
{
	double c13sq, sind, cosd, Jrr, Jvac;
	double Ue1sq, Ue2sq, Ue3sq, Um1sq, Um2sq, Um3sq, Ut1sq, Ut2sq, Ut3sq;
	double Lover4E, D21, D31;
	double sinD21, sinD31, sinD32;
	double sinsqD21_2, sinsqD31_2, sinsqD32_2, triple_sin;
	double Pme_TC, Pme_TV, Pmm, Pee;

	// ------------------------------------------------------------------------------------------------- //
	// First calculate useful simple functions of the oscillation parameters that don't depend on energy //
	// ------------------------------------------------------------------------------------------------- //
	c13sq = 1 - s13sq;

	// Ueisq's
	Ue3sq = s13sq;
	Ue2sq = c13sq * s12sq;

	// Umisq's, Utisq's and Jvac	 
	Um3sq = c13sq * s23sq;
	// Um2sq and Ut2sq are used here as temporary variables, will be properly defined later	 
	Ut2sq = s13sq * s12sq * s23sq;
	Um2sq = (1 - s12sq) * (1 - s23sq);
	  
	Jrr = sqrt(Um2sq * Ut2sq);
	sind = sin(delta);
	cosd = cos(delta);
	Um2sq = Um2sq + Ut2sq - 2 * Jrr * cosd;
	Jvac = 8 * Jrr * c13sq * sind;
	
	// ----------------------- //
	// Get all elements of Usq //
	// ----------------------- //
	Ue1sq = 1 - Ue3sq - Ue2sq;
	Um1sq = 1 - Um3sq - Um2sq;

	Ut3sq = 1 - Um3sq - Ue3sq;
	Ut2sq = 1 - Um2sq - Ue2sq;
	Ut1sq = 1 - Um1sq - Ue1sq;

	// --------------------------- //
	// Loop over the energy vector //
	// --------------------------- //
	for (unsigned int i = 0; i < E.size(); i++)
	{
		// ----------------------- //
		// Get the kinematic terms //
		// ----------------------- //
		Lover4E = eVsqkm_to_GeV_over4 * L / E[i];

		D21 = Dmsq21 * Lover4E;
		D31 = Dmsq31 * Lover4E;
		  
		sinD21 = sin(D21);
		sinD31 = sin(D31);
		sinD32 = sin(D31 - D21);

		triple_sin = sinD21 * sinD31 * sinD32;

		sinsqD21_2 = 2 * sinD21 * sinD21;
		sinsqD31_2 = 2 * sinD31 * sinD31;
		sinsqD32_2 = 2 * sinD32 * sinD32;

		// ----------------------------------------------------------------- //
		// Calculate the three necessary probabilities, separating TC and TV //
		// ----------------------------------------------------------------- //
		Pme_TC = (Ut3sq - Um2sq * Ue1sq - Um1sq * Ue2sq) * sinsqD21_2
			   + (Ut2sq - Um3sq * Ue1sq - Um1sq * Ue3sq) * sinsqD31_2
			   + (Ut1sq - Um3sq * Ue2sq - Um2sq * Ue3sq) * sinsqD32_2;
		
		Pme_TV = -Jvac * triple_sin;

		Pmm = 1 - 2 * (Um2sq * Um1sq * sinsqD21_2
					 + Um3sq * Um1sq * sinsqD31_2
					 + Um3sq * Um2sq * sinsqD32_2);

		Pee = 1 - 2 * (Ue2sq * Ue1sq * sinsqD21_2
					 + Ue3sq * Ue1sq * sinsqD31_2
					 + Ue3sq * Ue2sq * sinsqD32_2);

		// ---------------------------- //
		// Assign all the probabilities //
		// ---------------------------- //
		probs_returned[i][0][0] = Pee;														// Pee
		probs_returned[i][0][1] = Pme_TC - Pme_TV;											// Pem
		probs_returned[i][0][2] = 1 - Pee - probs_returned[i][0][1];  						// Pet

		probs_returned[i][1][0] = Pme_TC + Pme_TV;											// Pme
		probs_returned[i][1][1] = Pmm;														// Pmm
		probs_returned[i][1][2] = 1 - probs_returned[i][1][0] - Pmm;						// Pmt

		probs_returned[i][2][0] = 1 - Pee - probs_returned[i][1][0];						// Pte
		probs_returned[i][2][1] = 1 - probs_returned[i][0][1] - Pmm;						// Ptm
		probs_returned[i][2][2] = 1 - probs_returned[i][0][2] - probs_returned[i][1][2];	// Ptt
	} // i, energies
}
} // namespace NuFast

int main()
{
	int nE = 3; // number of energy points
	double L, Emin, Emax, Estep, rhoYe;
	vector<double> E(nE);
	vector<array<array<double, 3>, 3>> probs_returned(nE);
	double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31;
	int N_Newton;

	// ------------------------------- //
	// Set the experimental parameters //
	// ------------------------------- //
	L = 1300; // km
	Emin = 2; // GeV
	Emax = 5; // GeV
	Estep = (Emax - Emin) / nE; // GeV
	for (int i = 0; i < nE; i++)
		E[i] = Emin + i * Estep;
	rhoYe = 3 * 0.5; // g/cc

	// --------------------------------------------------------------------- //
	// Set the number of Newton-Raphson iterations which sets the precision. //
	// 0 is close to the single precision limit and is better than DUNE/HK   //
	// in the high statistics regime. Increasing N_Newton to 1,2,... rapidly //
	// improves the precision at a modest computational cost. Any negative   //
	// value (e.g. -1) accesses the (slow) exact analytic expression         //
	// --------------------------------------------------------------------- //
	N_Newton = 0;

	// ------------------------------------- //
	// Set the vacuum oscillation parameters //
	// ------------------------------------- //
	s12sq = 0.31;
	s13sq = 0.02;
	s23sq = 0.55;
	delta = 0.7 * M_PI;
	Dmsq21 = 7.5e-5; // eV^2
	Dmsq31 = 2.5e-3; // eV^2

	// ----------------------------------------------------------- //
	// Calculate all 9 oscillations probabilities for all energies //
	// ----------------------------------------------------------- //
	NuFast::Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, N_Newton, probs_returned);

	// --------------------------- //
	// Print out the probabilities //
	// --------------------------- //
	printf("L = %g (km), rhoYe = %g (g/cc)\n", L, rhoYe);
	printf("Probabilities:\n");
	printf("alpha beta E P(nu_alpha -> nu_beta)\n");
	for (int i = 0; i < nE; i++)
	{
		for (int alpha = 0; alpha < 3; alpha++)
		{
			for (int beta = 0; beta < 3; beta++)
			{
				printf("%d %d %g %.12f\n", alpha, beta, E[i], probs_returned[i][alpha][beta]);
			} // beta, 3
		} // alpha, 3
	} // i, nE, energy
}

