#include <cmath>
#include <cstdio>
#include <chrono>
#include <vector>
#include <array>

#include "Benchmark.h"
#include "NuFast_LBL.h"
#include "Exact_Cubic.h"
#include "Page.h"
#include "DMP.h"

using std::vector, std::array;

// Some constants
double const eVsqkm_to_GeV_over4 = 1e-9 / 1.97327e-7 * 1e3 / 4;
double const YerhoE2a = 1.52588e-4;
double const oneThird = 1. / 3;

// Compute the speed for nE energies, n times, and then calculate the mean over those n*nE probabilities
// N_Newton is used for NuFast only, ignored otherwise
double Speed_Helper(void Probability_Calculator(double, double, double, double, double, double, double, const vector<double> &, double, int, vector<array<array<double, 3>, 3>> &probs_returned), int n, int nE, int N_Newton)
{
	double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, Emin, Emax, Estep, rhoYe;
	vector<double> E(nE);
	vector<array<array<double, 3>, 3>> probs_returned(nE);

	// ------------------------------------- //
	// Set the vacuum oscillation parameters //
	// ------------------------------------- //
	s12sq = 0.31;
	s13sq = 0.02;
	s23sq = 0.55;
	delta = -0.7 * M_PI;
	Dmsq21 = 7.5e-5; // eV^2
	Dmsq31 = 2.5e-3; // eV^2

	// ------------------------------- //
	// Set the experimental parameters //
	// ------------------------------- //
	L = 1300; // km
	Emin = 0.5; // GeV
	Emax = 5;
	Estep = (Emax - Emin) / n;
	for (int i = 0; i < nE; i++)
		E[i] = Emin + i * Estep;
	rhoYe = 3 * 0.5; // g/cc

	std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < n; i++)
	{
		Probability_Calculator(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, N_Newton, probs_returned);
	} // i, E, n
	std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count() / (n * nE);
}

// Compute the speed many times and calculate the mean and standard deviation
// m is outer loop, n, is inner loop, nE is number of energy points
// N_Newton is used for NuFast only, ignored otherwise
void Speed(void Probability_Calculator(double, double, double, double, double, double, double, const vector<double> &, double, int, vector<array<array<double, 3>, 3>> &probs_returned), int m, int n, int nE, int N_Newton, double &mean, double &std)
{
	double s, speed_sum, speedsq_sum;

	speed_sum = 0;
	speedsq_sum = 0;
	for (int i = 0; i < m; i++)
	{
		s = Speed_Helper(Probability_Calculator, n, nE, N_Newton) * 1e9; // in ns
		speed_sum += s;
		speedsq_sum += s * s;
	} // i, m
	mean = speed_sum / m;
	std = sqrt(speedsq_sum / m - mean * mean);
}

int main()
{
	double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31;
	double L, rhoYe;
	vector<double> E;
	double Emin, Emax, Estep;
	double mean, std;
	int N_Newton, m, n, nE;
	vector<array<array<double, 3>, 3>> probs_returned_NuFast[3], probs_returned_Exact_Cubic, probs_returned_DMP;

	// ------------------------------------- //
	// Set the vacuum oscillation parameters //
	// ------------------------------------- //
	s12sq = 0.31;
	s13sq = 0.02;
	s23sq = 0.55;
	delta = -0.7 * M_PI;
	Dmsq21 = 7.5e-5;  // eV^2
	Dmsq31 = +2.5e-3; // eV^2

	// ------------------------------- //
	// Set the experimental parameters //
	// ------------------------------- //
	L = 1300; // km
	E = vector<double>(1, 2.5); // GeV
	rhoYe = 3 * 0.5; // g/cc

	N_Newton = 0; // <-- Sets the precision. 0 is close to the single precision limit and is better than DUNE/HK in the high statistics regime. Increasig N_Newton improves the precision at a modest computational cost

	// ------------------------------------------ //
	// Calculate all 9 oscillations probabilities //
	// ------------------------------------------ //
//	probs_returned_NuFast[0].resize(E.size());
	for (unsigned int i = 0; i < E.size(); i++)
		probs_returned_NuFast[0].emplace_back();
	NuFast::Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, N_Newton, probs_returned_NuFast[0]);
	probs_returned_Exact_Cubic.resize(E.size());
	Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, 0, probs_returned_Exact_Cubic);

	// ------------------------------- //
	// Print out the fractional errors //
	// ------------------------------- //
	printf("L = %g E = %g rhoYe = %g N_Newton = %d\n", L, E[0], rhoYe, N_Newton);
	printf("Precision (N_Newton = %d):\n", N_Newton);
	printf("alpha beta DeltaP/P\n");
	for (int alpha = 0; alpha < 3; alpha++)
	{
		for (int beta = 0; beta < 3; beta++)
		{
			printf("%d %d %g\n", alpha, beta, fabs(probs_returned_NuFast[0][0][alpha][beta] - probs_returned_Exact_Cubic[0][alpha][beta]) / probs_returned_Exact_Cubic[0][alpha][beta]);
		} // beta, 3
	} // alpha, 3

	// ----------------------------------- //
	// Write the fractional errors to file //
	// ----------------------------------- //
	nE = 1e3;
	// DUNE //
	Emin = 0.5;
	Emax = 5;
	Estep = (Emax - Emin) / nE;
	E.resize(nE + 1);
	for (int i = 0; i <= nE; i++)
		E[i] = Emin + i * Estep;
	probs_returned_Exact_Cubic.resize(E.size());
	Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, 0, probs_returned_Exact_Cubic);
	probs_returned_DMP.resize(E.size());
	Probability_Matter_LBL_DMP(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, 0, probs_returned_DMP);
	for (int j = 0; j < 2; j++)
	{
		probs_returned_NuFast[j].resize(E.size());
		NuFast::Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, j, probs_returned_NuFast[j]);
	} // j, N_Newton, 2
	FILE *dataf = fopen("data/Precision_DUNE.txt", "w");
	fprintf(dataf, "%g %g\n", L, rhoYe);
	for (int i = 0; i <= nE; i++)
	{
		fprintf(dataf, "%g ", E[i]);
		for (int j = 0; j < 2; j++)
		{
			fprintf(dataf, "%g %g ", fabs(probs_returned_NuFast[j][i][1][1] - probs_returned_Exact_Cubic[i][1][1]) / probs_returned_Exact_Cubic[i][1][1], fabs(probs_returned_NuFast[j][i][1][0] - probs_returned_Exact_Cubic[i][1][0]) / probs_returned_Exact_Cubic[i][1][0]);
		} // j, N_Newton, 2
		// now do DMP
		fprintf(dataf, "%g %g ", fabs(probs_returned_DMP[i][1][1] - probs_returned_Exact_Cubic[i][1][1]) / probs_returned_Exact_Cubic[i][1][1], fabs(probs_returned_DMP[i][1][0] - probs_returned_Exact_Cubic[i][1][0]) / probs_returned_Exact_Cubic[i][1][0]);
		fprintf(dataf, "\n");
	} // i, E, nE
	fclose(dataf);

	// HK //
	L = 295;
	Emin = 0.1;
	Emax = 2;
	Estep = (Emax - Emin) / nE;
	for (int i = 0; i <= nE; i++)
		E[i] = Emin + i * Estep;
	Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, 0, probs_returned_Exact_Cubic);
	for (int j = 0; j < 2; j++) // j, N_Newton
		NuFast::Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, j, probs_returned_NuFast[j]);
	dataf = fopen("data/Precision_HK.txt", "w");
	fprintf(dataf, "%g %g\n", L, rhoYe);
	for (int i = 0; i <= nE; i++)
	{
		fprintf(dataf, "%g ", E[i]);
		for (int j = 0; j < 2; j++)
		{
			fprintf(dataf, "%g %g ", fabs(probs_returned_NuFast[j][i][1][1] - probs_returned_Exact_Cubic[i][1][1]) / probs_returned_Exact_Cubic[i][1][1], fabs(probs_returned_NuFast[j][i][1][0] - probs_returned_Exact_Cubic[i][1][0]) / probs_returned_Exact_Cubic[i][1][0]);
		} // j, N_Newton, 2
		// now do DMP
		fprintf(dataf, "%g %g ", fabs(probs_returned_DMP[i][1][1] - probs_returned_Exact_Cubic[i][1][1]) / probs_returned_Exact_Cubic[i][1][1], fabs(probs_returned_DMP[i][1][0] - probs_returned_Exact_Cubic[i][1][0]) / probs_returned_Exact_Cubic[i][1][0]);
		fprintf(dataf, "\n");
	} // i, E, nE
	fclose(dataf);

	// JUNO //
	L = 50;
	rhoYe = 2.6 * 0.5;
	Emin = 1e-3;
	Emax = 10e-3;
	Estep = (Emax - Emin) / nE;
	for (int i = 0; i <= nE; i++)
		E[i] = Emin + i * Estep;
	Probability_Matter_LBL_Exact_Cubic(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, 0, probs_returned_Exact_Cubic);
	for (int j = 0; j < 2; j++) // j, N_Newton
		NuFast::Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rhoYe, j, probs_returned_NuFast[j]);
	dataf = fopen("data/Precision_JUNO.txt", "w");
	fprintf(dataf, "%g %g %d\n", L, rhoYe, N_Newton);
	for (int i = 0; i <= nE; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			fprintf(dataf, "%g %g\n", E[i], fabs(probs_returned_NuFast[j][i][0][0] - probs_returned_Exact_Cubic[i][0][0]) / probs_returned_Exact_Cubic[i][0][0]);
		} // j, N_Newton, 2
	} // i, E, nE
	fclose(dataf);

	// ---------- //
	// Speed test //
	// ---------- //
	dataf = fopen("data/Speed.txt", "w");
	m = int(1e2);
	printf("\nSpeed (this takes a moment):\n");
	for (int i = 0; i < 4; i++)
	{
		nE = pow(10, i);
		n = int(1e5 / nE);

		printf("Number of energy points: %d\n", nE);
		fprintf(dataf, "Number of energy points: %d\n", nE);

		// Vacuum probability
		printf("Vacuum (1/7)\n");
		Speed(NuFast::Probability_Vacuum_LBL, m, n, nE, 0, mean, std);
		printf("%g +- %g ns\n", mean, std);
		fprintf(dataf, "Vacuum %g += %g ns\n", mean, std);

		// NuFast at various N_Newton levels
		for (N_Newton = -1; N_Newton <= 3; N_Newton++)
		{
			printf("N_Newton = %d (%d/7)\n", N_Newton, N_Newton + 2);
			Speed(NuFast::Probability_Matter_LBL, m, n, nE, N_Newton, mean, std);
			printf("%g +- %g ns\n", mean, std);
			fprintf(dataf, "N_Newton=%d %g += %g ns\n", N_Newton, mean, std);
		} // N_Newton, 0, 3

		// Exact with cubic from Cardano/ZS
		printf("Exact_Cubic (6/7)\n");
		Speed(Probability_Matter_LBL_Exact_Cubic, m, n, nE, 0, mean, std);
		printf("%g +- %g ns\n", mean, std);
		fprintf(dataf, "Exact_Cubic %g += %g ns\n", mean, std);

		// Page's algorithm from 2309.06900
		printf("Page (7/7)\n");
		Speed(Probability_Matter_LBL_Page, m, n, nE, 0, mean, std);
		printf("%g +- %g ns\n\n", mean, std);
		fprintf(dataf, "Page %g += %g ns\n\n", mean, std);
	} // i, nE

	fclose(dataf);

	return 0;
}
