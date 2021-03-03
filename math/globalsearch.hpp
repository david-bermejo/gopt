#pragma once

#include "vector.hpp"
#include "hcint.hpp"

namespace gopt
{
	template <typename T, typename F, unsigned int N>
	Vector_t<T, N> globalsearch(F& f, const Vector_t<T, N>& lb, const Vector_t<T, N>& ub)
	{
		const T k = 25;

		const auto kexp = [&](const Vector_t<T, N>& x)
		{
			return std::exp(-k * f(x));
		};

		const T kexp_int = std::get<0>(hcint<1>(kexp, lb, ub));
		const auto minima_dist = [&](const Vector_t<T, N>& x)
		{
			return kexp(x) / kexp_int;
		};

		const auto global_finder = [&](const Vector_t<T, N>& x)
		{
			return x * minima_dist(x);
		};

		return std::get<1>(hcint<N>(global_finder, lb, ub));
	}
}