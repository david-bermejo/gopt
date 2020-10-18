#pragma once

#include <cmath>

namespace gopt
{
	template <typename T = double>
	constexpr T e = 2.7182818284590452353602875;

	template <typename T = double>
	constexpr T pi = 3.141592653589793238462643383279;

	template <typename T = double>
	constexpr T epsilon(const T x) { return std::nextafter(x, x + 0.1) - x; }
}