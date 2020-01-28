#pragma once

#include <cmath>

namespace math
{
	template <typename T>
	constexpr T epsilon = std::numeric_limits<T>::epsilon();
}