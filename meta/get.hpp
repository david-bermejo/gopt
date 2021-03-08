#pragma once

namespace meta
{
    struct no_type{};

    template <int N, typename... Args>
    struct get;

    template <int N, typename T, typename... Args>
    struct get<N, T, Args...>
    {
        using type = typename get<N-1, Args...>::type;
    };

    template <typename T, typename... Args>
    struct get<0, T, Args...>
    {
        using type = T;
    };
}