#pragma once

#include <utility>

namespace meta
{
	template <unsigned int N, unsigned int... Seq>
	constexpr std::integer_sequence<unsigned int, N+Seq...> add(std::integer_sequence<unsigned int, Seq...>) { return {}; }

	template<unsigned int Min, unsigned int Max>
	using make_integer_range = decltype(add<Min>(std::make_integer_sequence<unsigned int, Max-Min>()));
}