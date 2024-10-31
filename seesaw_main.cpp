//include files
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include "seesaw.h"

using namespace std;

int main()
{
	double p[Np] = {1.0, 2.0, 3.0, 			// mNu_light, dm21, dm3l (eV, eV^2, eV^2)
					4.0, 5.0, 6.0, 			// theta23, theta12, theta13
					0.1, 0.2, 0.3, 			// deltaCP, alpha1, alpha2
					0.4, 0.5, 0.6, 			// a1, a2, a3
					0.01, 0.02, 0.03, 		// phi1, phi2, phi3
					0.04, 0.05, 0.06}; 		// M1, M2, M3 (GeV)
	double obs[No] = {0.2, 0.4, 0.6};
	double err[No] = {0.1, 0.2, 0.3};
	double prior_mean[Np] = {1.0, -1.0, 0.0};
	double prior_std[Np] = {2.0, 0.5, 1.0};

	Seesaw model(p, obs, err, prior_mean, prior_std);
	double a;
	a = model.logposterior();
	cout << a << endl;

	cout << model.md21() << endl;
	cout << model.get_ordering() << endl;
	cout << model.Hmatrix12_real() << endl;

	double ps[Np] = {3., 4., -5.};
	model.set_params(ps);
	a = model.logposterior();
	cout << a << endl;

	cout << model.md21() << endl;
	cout << model.get_ordering() << endl;
	

	return(0);
}
