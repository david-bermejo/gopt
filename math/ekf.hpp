#pragma once

#include "algorithms.hpp"
#include "matrix.hpp"
#include "solvers.hpp"
#include "vector.hpp"

namespace gopt
{
	template <typename T, unsigned int Nx, unsigned int Nz, typename F, typename H>
	void ekf(Vector_t<T, Nx>& x,
			 Vector_t<T, Nz>& z,
			 Matrix_t<T, Nx, Nx>& P,
			 const Matrix_t<T, Nx, Nx>& Q,
			 const Matrix_t<T, Nz, Nz>& R,
			 F& f,
			 H& h,
			 const T& t0,
			 const T& tf)
	{
		const Vector_t<T, Nx> x_est = DP45(f, x, t0, tf, 1e-3, 1e-6, tf-t0);

		// Delta value to compute the Jacobians.
		const T dx = 1e-9;

		// Build Jacobian matrix of f and h.
		Matrix_t<T, Nx, Nx> Fk;
		Vector_t<T, Nx> f_eval[Nx];

		// Populate eval vector array.
		for (int i = 0; i < Nx; i++)
		{
			Vector_t<T, Nx> x_dx = x;
			x_dx[i] += dx;

			f_eval[i] = DP45(f, x_dx, t0, tf, 1e-3, 1e-6, tf-t0);
		}

		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Nx; j++)
				Fk[i][j] = (f_eval[j][i] - x_est[i]) / dx;

		// Build Jacobian matrix of f and h.
		Matrix_t<T, Nz, Nx> Hk;
		Vector_t<T, Nz> h_eval[Nx];
		const auto z_est = h(x_est, tf);

		// Populate eval vector array.
		for (int i = 0; i < Nx; i++)
		{
			Vector_t<T, Nx> x_dx = x_est;
			x_dx[i] += dx;

			h_eval[i] = h(x_dx, tf);
		}

		for (int i = 0; i < Nz; i++)
			for (int j = 0; j < Nx; j++)
				Hk[i][j] = (h_eval[j][i] - z_est[i]) / dx;

		const Matrix_t<T, Nx, Nx> P_est = Fk * P * transpose(Fk) + Q;
		const auto HT = transpose(Hk);
		const auto K = P_est * HT * inverse(Hk * P_est * HT + R);

		x = x_est + K * (z - z_est);
		P = (eye<T, Nx>() - K * Hk) * P_est;
	}
}