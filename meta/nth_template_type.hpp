#pragma once

#include "get.hpp"

namespace meta
{
    template <unsigned int N, typename T>
    struct nth_template_type;

    template <unsigned int N, template<typename...> class C, typename... Args>
    struct nth_template_type<N, C<Args...>>
    {
        using type = typename get<N-1, Args...>::type;
    };

    template <unsigned int N, typename T>
    using nth_template_type_t = typename nth_template_type<N, T>::type;
}