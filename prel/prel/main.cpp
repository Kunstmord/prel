//
//  main.cpp
//  prel
//
//  Created by George Oblapenko on 11/11/14.
//  Copyright (c) 2014 George Oblapenko. All rights reserved.
//

#include <iostream>
#include <kineticlib.h>

int main(int argc, const char * argv[]) {
	double start_T = 2000.; // the start temperature
	double end_T = 10000.0; // the end temperature
	double p = 100000.0; // atmospheric pressure
	double xN = 0.5; // the relative numeric density of atomic nitrogen
	double xN2 = 1. - xN;
	
	std::string cs_model = "VSS";  // the dissociation cross-section model
	bool vl_dependent = false;  // whether the dissociation cross-section is dependent on the vibrational level of the dissociating molecule

	int points_amt = 50; // the amount of points
	arma::vec T_array = arma::linspace<arma::vec>(start_T, end_T, points_amt); // start value of T, end value of T, amount of steps
	arma::vec prel = arma::zeros(50); // array of relaxation pressure values

	klib::MoleculeOneT N2 = klib::MoleculeOneT("N2");
	klib::Atom N = klib::Atom("N");
	arma::vec idata_N2N2 = klib::load_elastic_parameters(N2, N2); // N2 + N2 elastic interaction data
	arma::vec idata_N2N = klib::load_elastic_parameters(N2, N); // N2 + N elastic interaction data
	arma::vec ddata_N2N2 = klib::load_dissociation_parameters(N2, N2); // N2 + N2 dissociation data
	arma::vec ddata_N2N = klib::load_dissociation_parameters(N2, N); // N2 + N dissociation data
	arma::mat33 beta_matrix; // the matrix of beta integral brackets
	arma::vec3 right_parts; // the right-hand side of the system
	arma::vec3 results; // the values of the expansion coefficients g_{c,pq} - the order is g_{N2,10}, g_{N2,01}, g_{N,10}
	 
	double R_react_N2, R_react_N; // relaxation terms
	double T, n, rho, omega11_N2N, tau_rot_N2N, tau_rot_N2N2, eta_zeta_N2N2, eta_zeta_N2N, cV, SJsl, PJsl, SJslN2, PJslN2; // the numeric density of the mixture, the density of the mixture, Omega^{(1,1)}_{N2,N}, rotational relaxation times for N2+N and N2+N2, the quantity 4T \pi / (\eta_{cd} * \xi_{cd}), the quantity dE/dT
	double tmp_int00 = 0.0;
	beta_matrix.at(0, 0) = 1.5 * (1. - xN);
	beta_matrix.at(0, 2) = 1.5 * xN; // these appear due to the constraint conditions and are independent of temperature
	beta_matrix.zeros();
	right_parts[0] = 0.;
	for (int i = 0; i < points_amt; i++) {
		T = T_array[i];
		n = p / (KLIB_CONST_K * T);
		N2.renorm(T, xN2 * n);
		rho = xN2 * n * N2.mass + xN * n * N.mass;
		omega11_N2N = klib::omega(T, 1, 1, idata_N2N, "ESA", true, true); // we set nokt to true so that the result is of order 10e-07 instead of 10e-16
		tau_rot_N2N = klib::rot_rel_time_vss(T, idata_N2N, N2, n);
		tau_rot_N2N2 = klib::rot_rel_time_vss(T, idata_N2N2, N2, n);
		eta_zeta_N2N2 = 1. / (n * KLIB_CONST_K * tau_rot_N2N2);
		eta_zeta_N2N = 1. / (n * KLIB_CONST_K * tau_rot_N2N);
		cV = 1.5 * KLIB_CONST_K * n / rho + xN2 * (N2.mass * n / rho) * (N2.crot + N2.E_vibr_dT(T));
		beta_matrix.at(0, 1) = xN2 * klib::wt_poly_norm(T, N2);

		beta_matrix.at(1, 0) = -xN2 * N2.mass * N2.crot * (xN * (0.33333) * eta_zeta_N2N + xN2 * eta_zeta_N2N2); // does vibrational relaxation time play any significant role?
		beta_matrix.at(1, 1) = xN2 * N2.mass * N2.crot * (xN * eta_zeta_N2N + xN2 * eta_zeta_N2N2);
		beta_matrix.at(2, 0) = (2. / 9.) * xN2 * xN * (-16 * omega11_N2N + sqrt(1. / (KLIB_CONST_K * T)) * N2.mass * N2.crot * eta_zeta_N2N);
		beta_matrix.at(2, 2) = (2. / 9.) * xN2 * xN * (16 * omega11_N2N + sqrt(1. / (KLIB_CONST_K * T)) * N2.mass * N2.crot * eta_zeta_N2N);

		beta_matrix.at(1, 0) /= sqrt(KLIB_CONST_K * T);
		beta_matrix.at(1, 1) /= sqrt(KLIB_CONST_K * T);
		R_react_N2 = -klib::diss_rate_treanor_marrone(T, ddata_N2N, N2) * (xN * n) * (xN2 * n);
		
		R_react_N2 /= sqrt(KLIB_CONST_K * T);
		R_react_N = -2 * R_react_N2; // this is a binary mixture and the relaxation terms are related in a simple way

		SJsl = 0.0;
		PJsl = 0.0;
		SJslN2 = 0.0;
		PJslN2 = 0.0;

		for (int vl = 0; vl <= N2.num_vibr; vl++) {
			tmp_int00 = klib::diss_integral(T, 0, idata_N2N, N2, vl, true, vl_dependent, cs_model, true);
			SJsl += (1. / 3.) * (12 * tmp_int00  - klib::diss_integral(T, 1, idata_N2N, N2, vl, true, true, cs_model, true)) * N2.vibr_exp(T, vl) / N2.Z_vibr(T);
			PJsl += (N2.avg_vibr_energy(T, false) - N2.vibr[vl] / (KLIB_CONST_K * T)) * 8 * tmp_int00 * N2.vibr_exp(T, vl) / N2.Z_vibr(T);

			tmp_int00 = klib::diss_integral(T, 0, idata_N2N2, N2, vl, true, vl_dependent, cs_model, true);
			SJslN2 += (1. / 2.) * (12 * tmp_int00 - 8 * klib::diss_integral(T, 1, idata_N2N2, N2, vl, true, vl_dependent, cs_model, true)) * N2.vibr_exp(T, vl) / N2.Z_vibr(T);
			PJslN2 += (N2.avg_vibr_energy(T, false) - N2.vibr[vl] / (KLIB_CONST_K * T)) * 8 * tmp_int00 * N2.vibr_exp(T, vl) / N2.Z_vibr(T);
		}

		right_parts[1] = xN2 * (R_react_N2 * (1.5 * KLIB_CONST_K * T + N2.avg_full_energy(T) + N2.form)
			+ R_react_N * (1.5 * KLIB_CONST_K * T + N.form)) * klib::wt_poly_norm(T, N2) / (rho * T * cV) + PJsl * (xN * xN2) * n;
		right_parts[2] = xN * (R_react_N2 * (1.5 * KLIB_CONST_K * T + N2.avg_full_energy(T) + N2.form)
			+ R_react_N * (1.5 * KLIB_CONST_K * T + N.form)) * 1.5 / (rho * T * cV) - 2 * 1.5 * SJsl * (xN * xN2) * n;
		results = arma::solve(beta_matrix, right_parts);
		std::cout << "\nprel_{N2+N}=" << -KLIB_CONST_K * T * (xN2 * results[0] + xN * results[2]) * klib::Gamma_diss(T, N2, N, N, xN2 * n, xN * n, xN * n);

		R_react_N2 = klib::diss_rate_treanor_marrone(T, ddata_N2N2, N2) * (xN2 * n) * (xN2 * n);
		R_react_N2 /= sqrt(KLIB_CONST_K * T);
		R_react_N = -2 * R_react_N2;
		right_parts[1] = xN2 * (R_react_N2 * (1.5 * KLIB_CONST_K * T + N2.avg_full_energy(T) + N2.form)
			+ R_react_N * (1.5 * KLIB_CONST_K * T + N.form)) * klib::wt_poly_norm(T, N2) / (rho * T * cV) + PJslN2 * (xN2 * xN2) * n;
		right_parts[2] = xN * (R_react_N2 * (1.5 * KLIB_CONST_K * T + N2.avg_full_energy(T) + N2.form)
			+ R_react_N * (1.5 * KLIB_CONST_K * T + N.form)) * 1.5 / (rho * T * cV) - 2 * 1.5 * SJslN2 * (xN2 * xN2) * n;
		results = arma::solve(beta_matrix, right_parts);
		std::cout << "; prel_{N2+N2}=" << -KLIB_CONST_K * T * (xN2 * results[0] + xN * results[2]) * klib::Gamma_diss(T, N2, N, N, xN2 * n, xN * n, xN * n);
	}
}
