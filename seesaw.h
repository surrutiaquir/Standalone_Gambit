#ifndef SEESAW_H
#define SEESAW_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#define Np 18  // Number of model parameters
#define No 3  // NUmber of observables

#define Nmv 3 // Number of active (light) neutrino masses

class Seesaw
{
	private:
		// Model paramerters
    	double params[Np];
		// Experimental values of observables
    	double obs_ex[No];
		// Experimental uncertainty of observables
    	double err_ex[No];
		// Mean of Gaussian prior
		double prior_mean[Np];
		// Standard deviation of Gaussian prior
		double prior_std[Np];

		// Active (light) neutrino masses
		Eigen::Matrix3cd m_nu;
		double mv1, mv2, mv3;
		bool ordering;
		// Sterile (heavy) neutrino masses
		Eigen::Matrix3cd M_h;
		double M1, M2, M3;
		// PMNS matrix
		Eigen::Matrix3cd U_nu;
		double theta12, theta13, theta23, deltaCP, alpha1, alpha2;
		// R matrix
		Eigen::Matrix3cd R;
		// H, Theta & Vnu matrices
		Eigen::Matrix3cd H, Theta, Vnu;
		
	public:
		// Initializer
    	Seesaw(double p[Np], double obs[No], double err[No], double pr_mean[Np], double pr_std[Np]);

		// Set model parameters to new values
    	void set_params(double p[Np]);
		void set_masses();
		void set_PMNS();
		void set_Rmatrix();
		void set_H__Theta__Vnu();

		// Observable functions
    	double obs0();
		double obs1();
		double obs2();
		bool get_ordering();
		double md21();
		double get_mv1();
		double get_mv2();
		double get_mv3();
		double Hmatrix12_real();

		// Log likelihood
		double loglikelihood();

		// Log prior
		double logprior();

		// Log posterior
		double logposterior();

};

#endif
