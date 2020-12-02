#pragma once

#include "algorithms.hpp"
#include "constants.hpp"
#include "matrix.hpp"
#include "vector.hpp"

#include <array>

namespace gopt
{
	namespace internal
	{
		template <typename T, unsigned int N>
		std::array<Vector_t<T, N>, N + 1> build_simplex_vertices()
		{
			std::array<Vector_t<T, N>, N + 1> res;

			for (int i = 0; i < N + 1; i++)
			{
				res[i] = Vector_t<T, N>(0);
				for (int j = 0; j < i; j++)
					res[i][j] = -std::sqrt((N + 1) / static_cast<T>(N * (N - j + 2) * (N - j + 1)));
				if (i < N)
					res[i][i] = std::sqrt((N + 1) * (N - i + 1) / static_cast<T>(N * (N - i + 2)));
			}

			return res;
		}

		template <typename T>
		T poly_deriv_sqr(const T& t, unsigned int N)
		{
			const T beta = (N - 2) / static_cast<T>(2);
			const T tmp = 3 * t * t - 6 * (beta + 3) * t + 3 * (beta + 3) * (beta + 2);
			return tmp * tmp;
		}

		template <typename T>
		std::array<T, 3> poly_roots(unsigned int N)
		{
			std::array<T, 3> res;

			const T beta = (N - 2) / static_cast<T>(2);
			const T theta = std::acos(std::pow(beta + 3, static_cast<T>(-0.5)));

			res[0] = static_cast<T>(2) * std::sqrt(beta + 3) * std::cos((theta + 2 * pi<T>) / 3) + (beta + 3);
			res[1] = static_cast<T>(2) * std::sqrt(beta + 3) * std::cos((theta + 4 * pi<T>) / 3) + (beta + 3);
			res[2] = static_cast<T>(2) * std::sqrt(beta + 3) * std::cos(theta / 3) + (beta + 3);

			return res;
		}

		template <typename T>
		std::array<T, 3> poly_weights(unsigned int N)
		{
			std::array<T, 3> res;

			const std::array<T, 3> t = poly_roots<T>(N);
			const T beta = (N - 2) / static_cast<T>(2);

			for (int i = 0; i < 3; i++)
				res[i] = 6 * std::tgamma(beta + 4) / (2 * (N + 1) * std::tgamma(static_cast<T>(N) / 2) * t[i] * poly_deriv_sqr(t[i], N));

			return res;
		}
	}

	template <typename T, unsigned int Nx, unsigned int Nz, typename F, typename H>
	void ssrcqkf(Vector_t<T, Nx>& x,
				 Vector_t<T, Nz>& z,
				 Matrix_t<T, Nx, Nx>& Pxx,
				 const Matrix_t<T, Nx, Nx>& Q,
				 const Matrix_t<T, Nz, Nz>& R,
				 F& f,
				 H& h,
				 const T& t0,
				 const T& tf)
	{
		/*
		 * Calculate constants
		 */

		//N+1 vertices of the n-dimensional simplex.
		const static std::array<Vector_t<T, Nx>, Nx + 1> vertices = internal::build_simplex_vertices<T, Nx>();

		// 3rd degree Chevyshev-Laguerre polynomial roots.
		const static std::array<T, 3> poly_roots = internal::poly_roots<T>(Nx);

		// Sum weights.
		const static std::array<T, 3> weights = internal::poly_weights<T>(Nx);

		/*
		 * 1st step: State update
		 */

		// Calculate Cholesky factorization of Pxx.
		auto Cxx = cholesky(Pxx);

		// Calculate x_est using function f provided.
		Vector_t<T, Nx> x_est(0);
		Vector_t<T, Nx> X_est[3][2*Nx+2];

		for (int j = 0; j < 3; j++)
		{
			Vector_t<T, Nx> partial_v(0);

			for (int i = 0; i < Nx + 1; i++)
			{
				const static T roots[3] = { std::sqrt(2 * poly_roots[0]), std::sqrt(2 * poly_roots[1]), std::sqrt(2 * poly_roots[2]) };
				const Vector_t<T, Nx> tmp = roots[j] * (Cxx * vertices[i]);
				X_est[j][2 * i] = f(x + tmp);
				X_est[j][2 * i + 1] = f(x - tmp);
				partial_v += X_est[j][2 * i] + X_est[j][2 * i + 1];
			}

			x_est += partial_v * weights[j];
		}

		// Calculate Pxx_est using the computed vectors above.
		Matrix_t<T, Nx, Nx> Pxx_est(0);
		for (int j = 0; j < 3; j++)
		{
			Matrix_t<T, Nx, Nx> partial_m(0);
			for (int i = 0; i < 2 * Nx + 2; i++)
			{
				const Vector_t<T, Nx> tmp = X_est[j][i] - x_est;
				partial_m += outer(tmp, tmp);
			}

			Pxx_est += partial_m * weights[j];
		}

		// Add noise covariance matrix Q to Pxx_est
		Pxx_est += Q;

		/*
		 * 2nd step: Measurement update
		 */

		// Calculate Cholesky factorization of Pxx_est.
		std::cout << "Doing..." << std::endl;
		Cxx = cholesky(Pxx_est);
		std::cout << "Done." << std::endl;

		// Calculate z_est using function h provided.
		Vector_t<T, Nz> z_est(0);
		Vector_t<T, Nz> Z_est[3][2*Nx+2];
		//std::array<std::array<Vector_t<T, Nz>, 2 * Nx + 2>, 3> Z_est;

		for (int j = 0; j < 3; j++)
		{
			Vector_t<T, Nz> partial_v(0);

			for (int i = 0; i < Nx + 1; i++)
			{
				const static T roots[3] = { std::sqrt(2 * poly_roots[0]), std::sqrt(2 * poly_roots[1]), std::sqrt(2 * poly_roots[2]) };
				const Vector_t<T, Nx> tmp = roots[j] * (Cxx * vertices[i]);
				X_est[j][2 * i] = x + tmp;
				X_est[j][2 * i + 1] = x - tmp;
				Z_est[j][2 * i] = h(X_est[j][2 * i], tf);
				Z_est[j][2 * i + 1] = h(X_est[j][2 * i + 1], tf);
				partial_v += Z_est[j][2 * i] + Z_est[j][2 * i + 1];
			}

			z_est += partial_v * weights[j];
		}

		Matrix_t<T, Nx, Nz> Pxz(0);
		Matrix_t<T, Nz, Nz> Pzz(0);

		// Calculate Pxx_est using the computed vectors above.
		for (int j = 0; j < 3; j++)
		{
			Matrix_t<T, Nx, Nz> partial_m_xz(0);
			Matrix_t<T, Nz, Nz> partial_m_zz(0);

			for (int i = 0; i < 2 * Nx + 2; i++)
			{
				const Vector_t<T, Nz> tmp = Z_est[j][i] - z_est;
				partial_m_xz += outer(X_est[j][i] - x_est, tmp);
				partial_m_zz += outer(tmp, tmp);
			}

			Pxz += partial_m_xz * weights[j];
			Pzz += partial_m_zz * weights[j];
		}

		// Add noise covariance matrix R to Pzz
		Pzz += R;

		/*
		 * 3rd step: State Update
		 */

		// Pzz is always positive-definite, so Cholesky factorization can be used to invert it
		Matrix_t<T, Nx, Nz> K = Pxz * inverse(Pzz);
		x = x_est + K * (z - z_est);
		Pxx = Pxx_est - (K * Pzz * transpose(K));
	}
}