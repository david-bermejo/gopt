#pragma once

#include "algorithms.hpp"
#include "vector.hpp"
#include <variant>

namespace gopt
{
	// Ugly fix so that:
	// *  S=1: Ret = T;
	// * S!=1: Ret = Vector_t<T, N>;
	template <typename T, unsigned int S>
	using Ret = std::variant<T, Vector_t<T, S>>;

	// Integrate the function f over the given hypercube using a 3rd order Legendre cuadrature.
	template <unsigned int M = 1, typename T, typename F, unsigned int N>
	Ret<T, M> hcint(F& f, const Vector_t<T, N>& lb, const Vector_t<T, N>& ub)
	{
		const std::array<T, 3> w { (T)5/9, (T)8/9, (T)5/9 };
		const std::array<T, 3> eps { -std::sqrt((T)3/5), 0, std::sqrt((T)3/5) };

		Ret<T, M> sum;

		if constexpr (M == 1)
			sum = (T)0;
		else
			sum = Vector_t<T, M>(0);

		constexpr int S = w.size();
		const int size = std::pow(S, N);

		#pragma omp parallel for
		for (int i = 0; i < size; i++)
		{
			Vector_t<T, N> x_i;
			T weight = 1;
			int divider = 1;

			for (int j = N-1; j >= 0; j--)
			{
				const int index = (i / divider) % S;

				weight *= w[index];
				x_i[j] = (ub[j] - lb[j]) / 2 * eps[index] + (ub[j] + lb[j]) / 2;
				divider *= S;
			}

			const Ret<T, M> tmp = weight * f(x_i);

			#pragma omp critical
			{ sum = std::get<M != 1>(sum) + std::get<M != 1>(tmp); }
		}

		return sum;
	}
}