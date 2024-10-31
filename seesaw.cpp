//include files
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <complex>
#include <vector>
#include "seesaw.h"
#include "NeutrinoInterpolator.hpp"


Seesaw::Seesaw(double p[Np], double obs[No], double err[No], double pr_mean[Np], double pr_std[Np])
{
	int Nparams = Np;
	for(int i=0; i<Nparams; i++)
	{
		params[i] = p[i];
		prior_mean[i] = pr_mean[i];
		prior_std[i] = pr_std[i];
	}

	set_masses();
	set_PMNS();
	set_Rmatrix();
	set_H__Theta__Vnu();

	if (ordering)
	{
		NeutrinoInterpolator spline_t12_n("data/NuFit_v4.1/T12n.csv");
		std::cout << "NI: " << spline_t12_n.eval(0.2) << std::endl;
	}
	else
	{
		NeutrinoInterpolator spline_t12_i("data/NuFit_v4.1/T23i.csv");
		std::cout << "NI: " << spline_t12_i.eval(0.2) << std::endl;
	}

	int Nobs = No;
	for(int i=0; i<Nobs; i++)
	{
		obs_ex[i] = obs[i];
		err_ex[i] = err[i];
	}
}

void Seesaw::set_masses()
{
	// First three parameters are mNu_light (eV), dm21 (eV^2), dm3l (eV^2)
	double mNu_light = params[0];
	double dm21 = params[1];
	double dm3l = params[2];
	if (dm3l > 0.) // Normal ordering
	{
		ordering = true; //true for normal ordering, false for inverted ordering.
		mv1 = mNu_light * 1e-9; // Convert to GeV
		mv2 = pow(mNu_light * mNu_light + dm21, 0.5) * 1e-9; // Convert to GeV
		mv3 = pow(mNu_light * mNu_light + dm3l, 0.5) * 1e-9; // Convert to GeV
	}
	else // Inverted ordering
	{
		ordering = false; //true for normal ordering, false for inverted ordering.
		mv3 = mNu_light * 1e-9; // Convert to GeV
		mv2 = pow(mNu_light * mNu_light - dm3l, 0.5) * 1e-9; // Convert to GeV
		mv1 = pow(mNu_light * mNu_light - dm3l - dm21, 0.5) * 1e-9; // Convert to GeV
	}
	m_nu(0,1) = 0.0;
	m_nu(0,2) = 0.0;
	m_nu(1,0) = 0.0;
	m_nu(1,2) = 0.0;
	m_nu(2,0) = 0.0;
	m_nu(2,1) = 0.0;
	m_nu(0,0) = mv1;
	m_nu(1,1) = mv2;
	m_nu(2,2) = mv3;

	// Parameters 16-18 are M1, M2, M3 (GeV)
	M1 = params[15];
	M2 = params[16];
	M3 = params[17];
	M_h(0,1) = 0.0;
	M_h(0,2) = 0.0;
	M_h(1,0) = 0.0;
	M_h(1,2) = 0.0;
	M_h(2,0) = 0.0;
	M_h(2,1) = 0.0;
	M_h(0,0) = M1;
	M_h(1,1) = M2;
	M_h(2,2) = M3;
}

void Seesaw::set_PMNS()
{
	// Parameters 4-9 are theta23, theta12, theta13, delta, alpha1, alpha2
	theta23 = params[3];
	theta12 = params[4];
	theta13 = params[5];
	deltaCP = params[6];
	alpha1 = params[7];
	alpha2 = params[8];
	std::complex<double> I(0.0, 1.0);
	Eigen::Matrix3cd V_23, V_13, V_12, U_pd, U_nd, Maj_phase;
	V_23 << 1.0, 0.0, 0.0,
            0.0, cos(theta23), sin(theta23),
            0.0, -sin(theta23), cos(theta23);
	V_13 << cos(theta13), 0.0, sin(theta13),
			0.0, 1.0, 0.0,
			-sin(theta13), 0.0, cos(theta13);
	V_12 << cos(theta12), sin(theta12), 0.0,
			-sin(theta12), cos(theta12), 0.0,
			0.0, 0.0, 1.0;
	U_pd << exp(-I*deltaCP/2.0), 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, exp(I*deltaCP/2.0);
	U_nd << exp(I*deltaCP/2.0), 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, exp(-I*deltaCP/2.0);
	Maj_phase << 1.0, 0.0, 0.0,
				0.0, exp(I*alpha1), 0.0,
				0.0, 0.0, exp(I*alpha2);
	U_nu = V_23 * U_pd * V_13 * U_nd * V_12 * Maj_phase;
}

void Seesaw::set_Rmatrix()
{
	// Parameters 10-15 are a1, a2, a3, phi1, phi2, phi3
	double a1 = params[9];
	double a2 = params[10];
	double a3 = params[11];
	double phi1 = params[12];
	double phi2 = params[13];
	double phi3 = params[14];
	std::complex<double> I(0.0, 1.0);
	Eigen::Matrix3cd A, Phi, R_O, R_H;
	A(0,0) = 0.0;
	A(0,1) = a1;
	A(0,2) = a2;
	A(1,0) = -a1;
	A(1,1) = 0.0;
	A(1,2) = a3;
	A(2,0) = -a2;
	A(2,1) = -a3;
	A(2,2) = 0.0;
	Phi(0,0) = I*0.0;
	Phi(0,1) = I*phi1;
	Phi(0,2) = I*phi2;
	Phi(1,0) = -I*phi1;
	Phi(1,1) = I*0.0;
	Phi(1,2) = I*phi3;
	Phi(2,0) = -I*phi2;
	Phi(2,1) = -I*phi3;
	Phi(2,2) = I*0.0;
	R_O = A.exp();
    R_H = Phi.exp();
	R = R_O * R_H;
}

void Seesaw::set_H__Theta__Vnu()
{
	Eigen::Matrix3cd One = Eigen::Matrix3cd::Identity();
	Eigen::Matrix3cd sqrt_mnu; // sqrt(m_nu)
	Eigen::Matrix3cd H2inv;  // H^-2
	std::complex<double> I(0.0, 1.0);

	sqrt_mnu << sqrt(m_nu(0,0).real()), 0.0, 0.0,
                0.0, sqrt(m_nu(1,1).real()), 0.0,
                0.0, 0.0, sqrt(m_nu(2,2).real());

	H2inv = One + sqrt_mnu * R * M_h.inverse() * R.adjoint() * sqrt_mnu;
	H = H2inv.sqrt().inverse();

	Theta = I * U_nu * H * m_nu.sqrt() * R * M_h.sqrt().inverse(); // Uas
	Vnu = U_nu * H; // Uaa
}

void Seesaw::set_params(double p[Np])
{
	int Nparams = Np;
	for(int i=0; i<Nparams; i++)
	{
		params[i] = p[i];
	}
	set_masses();
	set_PMNS();
	set_Rmatrix();
	set_H__Theta__Vnu();
}

double Seesaw::obs0()
{
	double o = 0;
	o = params[0] * params[1];
	return o;
}

double Seesaw::obs1()
{
	double o = 0;
	o = sin(params[2]) * pow(params[0], 2);
	return o;
}

double Seesaw::obs2()
{
	double o = 0;
	o = tanh(params[1]) * params[0] - params[2];
	return o;
}

double Seesaw::loglikelihood()
{
	double l = 0;
	l += pow(obs0() - obs_ex[0], 2) / pow(err_ex[0], 2);
	l += pow(obs1() - obs_ex[1], 2) / pow(err_ex[1], 2);
	l += pow(obs2() - obs_ex[2], 2) / pow(err_ex[2], 2);

	return l;
}

double Seesaw::logprior()
{
	double p = 0;
	for(int i=1; i<2; i++)
	{
		p += pow(params[i] - prior_mean[i], 2) / pow(prior_std[i], 2);
	}

	return p;
}

double Seesaw::logposterior()
{
	double l = 0;
	double p = 0;
	l = loglikelihood(); 
	p = logprior();
	return l + p;
}

bool Seesaw::get_ordering()
{
    return (mv3 >= mv1); // Returns true for normal ordering, false for inverted ordering.
}

double Seesaw::md21()
{
    return pow(mv2, 2) - pow(mv1, 2);
}

double Seesaw::get_mv1()
{
	return mv1;
}

double Seesaw::get_mv2()
{
	return mv2;
}

double Seesaw::get_mv3()
{
	return mv3;
}

double Seesaw::Hmatrix12_real()
{
	return H(0, 1).real();
}


