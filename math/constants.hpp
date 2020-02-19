#pragma once

#include <cmath>

namespace math
{
	template <typename T>
	constexpr T e = 2.7182818284590452353602875;

	template <typename T>
	constexpr T pi = 3.141592653589793238462643383279;

	template <typename T>
	constexpr T epsilon = std::numeric_limits<T>::epsilon();
}