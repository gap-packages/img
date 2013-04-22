//  (C) Copyright Howard Hinnant 2005-2011.
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).
//
//  See http://www.boost.org/libs/type_traits for most recent version including documentation.
//  See for GPL compatibility  http://www.gnu.org/licenses/license-list.html 

//  Details are in namespace detail.  Every effort has been made to make
//  combine_discontinuous and permute as fast as possible.  They minimize the number
//  of swaps that are performed. Everything else builds on these two primitives. 
//  The most complicated algorithm is for_each_reversible_permutation.  But it
//  builds on combine_discontinuous and permute and I believe represents a minimum
//  number of swaps.  Without care, algorithms such as for_each_reversible_permutation
//  will take longer than for_each_permutation instead of the intended half the time.

//  Speed is everything.  Lest you could just use std::next_permutation and manually
//  eliminate duplicate permutations.  If the implementation fails in being orders
//  of magnitude faster than that, then it has failed miserably.

#pragma once


#include <iterator>
#include <algorithm>
#include <limits>
#include <stdexcept>

// from http://marknelson.us/2002/03/01/next-permutation/

template<typename TInt>
inline TInt naive_next_combination(TInt N, TInt R, TInt c[])
{
    TInt i = 0;
    while ( i <R-1 && c[i+1]==c[i]+1) // for each bump
        //c[i] = i++;                 // fall back
        c[i] = i;                 
        i++;
    return N - ++c[i];              // push forward and verify
}
 
/// todo: need to test because of changes.
template<typename TInt>
inline TInt naive_next_combination_vec(TInt N, TInt R, std::vector<TInt> & c)
{
    size_t size = c.size();
    size_t i = 0;
    while ( i <size-1 && c[i+1]==c[i]+1) // for each bump
    {
        c[i] = i;                 // fall back
        i++;
    }
    return N - ++c[i];              // push forward and verify
}
