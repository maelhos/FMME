#pragma once

#include <NTL/vector.h>
#include "lzz_pXi.h"

NTL_CLIENT

template <size_t N>
using NVec = zz_pXi<N>;

template <size_t N>
zz_p getEval(const NVec<N>& vec, const array<zz_p, N>& point)
{
    return _getEval_impl(vec, point.data());
}

template <size_t N>
inline zz_p _getEval_impl(const NVec<N>& vec, const zz_p* point)
{
    return _getEval_impl(vec.rep[point[N - 1]._zz_p__rep], point);
}

inline zz_p _getEval_impl(const NVec<1>& vec, const zz_p* point)
{
    return vec.rep[point[0]._zz_p__rep];
}