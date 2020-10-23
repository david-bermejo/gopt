#pragma once

#include <cmath>

namespace gopt
{
	template <typename T = double>
	constexpr T e = static_cast<T>(2.7182818284590452353602875);

	template <typename T = double>
	constexpr T pi = static_cast<T>(3.141592653589793238462643383279);

	template <typename T = double>
	constexpr inline T epsilon(const T& x) { return std::nextafter(x, x + 0.1) - x; }

	template <typename T = double>
	constexpr inline T weak_epsilon = static_cast<T>(1.0e-7);

	template <typename T = double>
	constexpr inline T radians(const T& deg) { return deg * pi<T> / static_cast<T>(180); }

	template <typename T = double>
	constexpr inline T degrees(const T& rad) { return rad * static_cast<T>(180) / pi<T>; }
}