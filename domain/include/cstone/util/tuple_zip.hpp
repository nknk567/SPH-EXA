//
// Created by Noah Kubli on 29.08.23.
//

#pragma once

#include <tuple>

template <typename Tuple1, typename Tuple2, size_t ...Is>
decltype(auto) zip_tuples_impl(const Tuple1& tuple1, const Tuple2 &tuple2, std::index_sequence<Is...> is)
{
    return std::make_tuple(std::tie(std::get<Is>(tuple1), std::get<Is>(tuple2))...);
}

template <typename Tuple1, typename Tuple2>
decltype(auto) zip_tuples(const Tuple1& tuple1, const Tuple2 &tuple2)
{
    auto seq = std::make_index_sequence<std::tuple_size_v<std::decay_t<Tuple1>>>();
    return zip_tuples_impl(tuple1, tuple2, seq);
}
