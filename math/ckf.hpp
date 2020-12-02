#pragma once

#include "algorithms.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "solvers.hpp"
#include <array>

namespace gopt
{
	template <typename T, unsigned int Nx, unsigned int Nz, typename F, typename H>
	void ckf(Vector_t<T, Nx>& x,
			 const Vector_t<T, Nz>& z,
			 Matrix_t<T, Nx, Nx>& P,
			 const Matrix_t<T, Nx, Nx>& Q,
			 const Matrix_t<T, Nz, Nz>& R,
			 F& f,
			 H& h,
			 const T& t0,
			 const T& tf)
	{
		static const std::array<Vector_t<T, Nx>, Nx> vertices = []() { std::array<Vector_t<T, Nx>, Nx> res; const T v = std::sqrt(static_cast<T>(Nx)); for (int i = 0; i < Nx; i++) { Vector_t<T, Nx> tmp(0); tmp[i] = v; res[i] = tmp; } return res; }();
		
		// Estimation update
		const auto C = cholesky(P);
		Vector_t<T, Nx> x_est(0);
		Matrix_t<T, Nx, Nx> P_est(0);

		for (int i = 0; i < Nx; i++)
		{
			const Vector_t<T, Nx> tmp = C * vertices[i];

			const Vector_t<T, Nx> x_eval_p = DP45(f, x + tmp, t0, tf, 1e-9, 1e-9);
			const Vector_t<T, Nx> x_eval_m = DP45(f, x - tmp, t0, tf, 1e-9, 1e-9);

			x_est += x_eval_p + x_eval_m;
			P_est += outer(x_eval_p, x_eval_p) + outer(x_eval_m, x_eval_m);
		}
		x_est /= (2 * Nx);
		P_est /= (2 * Nx);
		P_est += Q - outer(x_est, x_est);

		// Measurement update
		const auto C_est = cholesky(P_est);
		Vector_t<T, Nz> z_est(0);
		Matrix_t<T, Nx, Nz> Pxz(0);
		Matrix_t<T, Nz, Nz> Pzz(0);

		for (int i = 0; i < Nx; i++)
		{
			const Vector_t<T, Nx> tmp = C_est * vertices[i];

			const Vector_t<T, Nx> x_eval_p = x_est + tmp;
			const Vector_t<T, Nx> x_eval_m = x_est - tmp;

			const Vector_t<T, Nz> z_eval_p = h(x_eval_p, tf);
			const Vector_t<T, Nz> z_eval_m = h(x_eval_m, tf);

			z_est += z_eval_p + z_eval_m;
			Pzz += outer(z_eval_p, z_eval_p) + outer(z_eval_m, z_eval_m);
			Pxz += outer(x_eval_p, z_eval_p) + outer(x_eval_m, z_eval_m);
		}
		z_est /= (2 * Nx);
		Pzz /= (2 * Nx);
		Pzz += R - outer(z_est, z_est);
		Pxz /= (2 * Nx);
		Pxz -= outer(x_est, z_est);

		const auto K = Pxz * inverse(Pzz);
		x = x_est + K*(z - z_est);
		P = P_est - K * Pzz * transpose(K);
	}
}