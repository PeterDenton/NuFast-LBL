#include <cmath>
#include <vector>
#include <array>

#include "DMP.h"
#include "Benchmark.h"

using std::vector, std::array;

#define sq(x) ((x)*(x))

void Probability_Matter_LBL_DMP(double s12sq, double s13sq, double s23sq, double delta, double Dmsq21, double Dmsq31, double L, vector<double> E, double rhoYe, int empty, vector<array<array<double, 3>, 3>> *probs_returned)
{
	double Dmsqeea, c2phi, a12, cphisq, sphisq, s2phi, cphi13sq, Dl21, c2psi, cpsisq, spsisq, Dl31, L4E, Delta21, Delta31, Delta32, sDelta21, sDelta31, sDelta32, dab, Jrm, C31, C32, C21, D, tmp;

	double c12, s12, s212, c212, c12sq, c13, s13, s213, c213, c13sq, c23, s23, c23sq, Dmsqee, cd, sd, a, Jrrm;
	c12sq = 1 - s12sq;
	c212 = c12sq - s12sq;
	s12 = sqrt(s12sq);
	c12 = sqrt(c12sq);
	s212 = 2 * s12 * c12;
	c23sq = 1 - s23sq;
	c23 = sqrt(c23sq);
	s23 = sqrt(s23sq);
	c13sq = 1 - s13sq;
	c13 = sqrt(c13sq);
	s13 = sqrt(s13sq);
	s213 = 2 * s13 * c13;
	c213 = c13sq - s13sq;
	Dmsqee = Dmsq31 - s12sq * Dmsq21;
	cd = cos(delta);
	sd = sin(delta);

	for (unsigned int i = 0; i < E.size(); i++)
	{
		a = rhoYe * E[i] * YerhoE2a;

		tmp = c213 - a / Dmsqee;
		Dmsqeea = Dmsqee * sqrt(sq(tmp) + sq(s213));
		c2phi = (Dmsqee * c213 - a) / Dmsqeea;
		a12 = 0.5 * (a + Dmsqee - Dmsqeea);

		cphisq = 0.5 * (1 + c2phi);
		sphisq = 0.5 * (1 - c2phi);
		s2phi = s213 * Dmsqee / Dmsqeea;
		cphi13sq = cphisq * c13sq + sphisq * s13sq + s2phi * c13 * s13;

		tmp = c212 - a12 / Dmsq21;
		Dl21 = Dmsq21 * sqrt(sq(tmp) + cphi13sq * sq(s212));
		c2psi = (Dmsq21 * c212 - a12) / Dl21;
		cpsisq = 0.5 * (1 + c2psi);
		spsisq = 0.5 * (1 - c2psi);

		Dl31 = Dmsq31 + 0.25 * a + 0.5 * (Dl21 - Dmsq21) + 0.75 * (Dmsqeea - Dmsqee);

		Jrm = s23 * c23 * cphisq * sqrt(sphisq * cpsisq * spsisq);
		Jrrm = Jrm / cphisq;

		L4E = eVsqkm_to_GeV_over4 * L / E[i];
		Delta21 = Dl21 * L4E;
		Delta31 = Dl31 * L4E;
		Delta32 = Delta31 - Delta21;

		sDelta21 = sin(Delta21);
		sDelta31 = sin(Delta31);
		sDelta32 = sin(Delta32);

		// Pme
		dab = 0;
		C31 = s23sq * sphisq * cphisq * cpsisq + Jrm * cd;
		C32 = s23sq * sphisq * cphisq * spsisq - Jrm * cd;
		C21 = cphisq * spsisq * cpsisq * (c23sq - sphisq * s23sq) + Jrm * cd * c2psi;
		D = -Jrm * sd;
		(*probs_returned)[i][1][0] = dab + 4 * (C31 * sq(sDelta31) + C32 * sq(sDelta32) + C21 * sq(sDelta21)) + 8 * D * sDelta21 * sDelta31 * sDelta32;

		// Pmm
		dab = 1;
		C31 = -cphisq * s23sq * (c23sq * spsisq + s23sq * sphisq * cpsisq) - 2 * s23sq * Jrm * cd;
		C32 = -cphisq * s23sq * (c23sq * cpsisq + s23sq * sphisq * spsisq) + 2 * s23sq * Jrm * cd;
		C21 = -(c23sq * cpsisq + s23sq * sphisq * spsisq) * (c23sq * spsisq + s23sq * sphisq * cpsisq) - 2 * (c23sq - sphisq * s23sq) * c2psi * Jrrm * cd + sq(2 * Jrrm * cd);
		D = 0;
		(*probs_returned)[i][1][1] = dab + 4 * (C31 * sq(sDelta31) + C32 * sq(sDelta32) + C21 * sq(sDelta21)) + 8 * D * sDelta21 * sDelta31 * sDelta32; // me
	} // i, energies
}
