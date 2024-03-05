#pragma once
#include <iostream>
#include "Data.h"

class FiniteDifferenceExample {
	double rho;
	double u;
	double Gamma;
	double L;
	double Pe;
	double phi_0;
	double phi_L;
	int N;
	std::string mesh;
	double r_e;
	std::string convection;
	void getExactSolution(Data& data) {
		data.x_exact.resize(101);
		for (int i = 0; i < 101; i++) {
			data.x_exact[i] = i * (L / 100);
		}
		data.phi_exact.resize(101);
		for (int i = 0; i < 101; i++) {
			data.phi_exact[i] = phi_0 + (exp(data.x_exact[i] * Pe / L) - 1) / (exp(Pe) - 1) * (phi_L - phi_0);
		}
	}
	void generateMesh(Data& data) {
		data.x_calculated.resize(N + 1);
		if (mesh == "uniform") {
			for (int i = 0; i < N + 1; i++) {
				data.x_calculated[i] = i * (L / N);
			}
		}
		else if (mesh == "non-uniform") {
			double dx = L * (1 - r_e) / (1 - pow(r_e, N));
			data.x_calculated[N] = L;
			if (u > 0) {
				for (int i = 1; i < N; i++) {
					data.x_calculated[i] = data.x_calculated[i - 1] + dx;
					dx = r_e * dx;
				}
			}
			else if (u < 0) {
				for (int i = N - 1; 0 < i; i--) {
					data.x_calculated[i] = data.x_calculated[i + 1] - dx;
					dx = r_e * dx;
				}
			}
		}
	}
	void discretize(Data& data, std::vector<double>& A_W, std::vector<double>& A_P, std::vector<double>& A_E, std::vector<double>& b) {
		std::vector<double> A_W_d(N - 1), A_P_d(N - 1), A_E_d(N - 1), A_W_c(N - 1), A_P_c(N - 1), A_E_c(N - 1);
		for (int i = 0; i < N - 1; i++) {
			A_W_d[i] = -2 * Gamma / ((data.x_calculated[i + 2] - data.x_calculated[i]) * (data.x_calculated[i + 1] - data.x_calculated[i]));
		}
		for (int i = 0; i < N - 1; i++) {
			A_E_d[i] = -2 * Gamma / ((data.x_calculated[i + 2] - data.x_calculated[i]) * (data.x_calculated[i + 2] - data.x_calculated[i + 1]));
		}
		for (int i = 0; i < N - 1; i++) {
			A_P_d[i] = -(A_W_d[i] + A_E_d[i]);
		}
		if (convection == "UDS") {
			for (int i = 0; i < N - 1; i++) {
				A_W_c[i] = -std::max(rho * u, 0.0) / (data.x_calculated[i + 1] - data.x_calculated[i]);
			}
			for (int i = 0; i < N - 1; i++) {
				A_E_c[i] = std::min(rho * u, 0.0) / (data.x_calculated[i + 2] - data.x_calculated[i + 1]);
			}
			for (int i = 0; i < N - 1; i++) {
				A_P_c[i] = -(A_W_c[i] + A_E_c[i]);
			}
		}
		else if (convection == "CDS") {
			for (int i = 0; i < N - 1; i++) {
				A_W_c[i] = -rho * u / (data.x_calculated[i + 2] - data.x_calculated[i]);
			}
			for (int i = 0; i < N - 1; i++) {
				A_E_c[i] = rho * u / (data.x_calculated[i + 2] - data.x_calculated[i]);
			}
		}
		for (int i = 0; i < N - 1; i++) {
			A_W[i] = A_W_d[i] + A_W_c[i];
		}
		for (int i = 0; i < N - 1; i++) {
			A_P[i] = A_P_d[i] + A_P_c[i];
		}
		for (int i = 0; i < N - 1; i++) {
			A_E[i] = A_E_d[i] + A_E_c[i];
		}
		b[0] = -A_W[0] * phi_0;
		b[N - 2] = -A_E[N - 2] * phi_L;
	}
	void Thomas(std::vector<double>& L, std::vector<double>& D, std::vector<double>& U, std::vector<double>& b) {
		for (int i = 0; i < D.size() - 1; i++) {
			U[i] = U[i] / D[i];
			D[i + 1] = D[i + 1] - L[i + 1] * U[i];
		}
		b[0] = b[0] / D[0];
		for (int i = 1; i < D.size(); i++) {
			b[i] = (b[i] - L[i] * b[i - 1]) / D[i];
		}
		for (int i = D.size() - 2; 0 <= i; i--) {
			b[i] = b[i] - U[i] * b[i + 1];
		}
	}
	void save(Data& data, std::vector<double>& b) {
		data.phi_calculated.resize(N + 1);
		data.phi_calculated[0] = phi_0;
		for (int i = 1; i < N; i++) {
			data.phi_calculated[i] = b[i - 1];
		}
		data.phi_calculated[N] = phi_L;
	}
public:
	FiniteDifferenceExample(double rho, double u, double Gamma, double L, double phi_0, double phi_L, int N, std::string mesh, double r_e, std::string convection) :rho(rho), u(u), Gamma(Gamma), L(L), Pe(rho* u* L / Gamma), phi_0(phi_0), phi_L(phi_L), N(N), mesh(mesh), r_e(r_e), convection(convection) {}
	void solve(Data& data) {
		std::vector<double> A_W(N - 1), A_P(N - 1), A_E(N - 1), b(N - 1);
		getExactSolution(data);
		generateMesh(data);
		discretize(data, A_W, A_P, A_E, b);
		Thomas(A_W, A_P, A_E, b);
		save(data, b);
	}
};
