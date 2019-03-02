#ifndef ARGSORT_HPP_FF5D00F6
#define ARGSORT_HPP_FF5D00F6

#include <algorithm>
#include <iterator>

/*! \file
\brief C++ port of Numpy's `argsort`.

Provides something similar to Numpy's `argsort`, i.e. computing the indices that
sort an array without actually sorting the array, with the style of the C++
standard library. Example:
~~~{.cpp}
#include <vector>
#include "argsort.h"
int main() {
    std::vector<int> v {3, 1, 2};
    
    std::vector<size_t> sorted_indices(v.size());
    argsort::argsort(sorted_indices.begin(), sorted_indices.end(), v.begin());
    // now `sorted_indices` is {1, 2, 0}
    
    std::vector<int> sorted_v(v.size());
    argsort::apply(
        sorted_indices.begin(), sorted_indices.end(),
        v.begin(),
        sorted_v.begin()
    );
    // now `sorted_v` is {v[1], v[2], v[0]} = {1, 2, 3}
}
~~~
*/

/*!
\brief Namespace for all the functions in `argsort.hpp`.
*/
namespace argsort {
    /*!
    \brief Fill a sequence of length `n` with indices `0...(n-1)`.
    
    Fills the sequence [`begin`, `end`) with increasing values starting from the
    default constructed value of the referenced type and increasing with `++`.
    */
    template<typename Iterator>
    void fill_range(Iterator begin, Iterator end) {
        using Type = typename std::iterator_traits<Iterator>::value_type;
        for (Type i {}; begin != end; ++begin) {
            *begin = i;
            ++i;
        }
    }
    
    /*!
    \brief Compute indices that indicate how to sort an array.
    
    Fills the sequence [`indices_begin`, `indices_end`) with indices such that,
    if the indices are `{i0, i1, ..., iN}` then `{*(data + i0), *(data + i1),
    ..., *(data + iN)}` is sorted. The length of the sequence starting at `data`
    must be at least the same of [`indices_begin`, `indices_end`).
    */
    template<typename IndexIterator, typename ValueIterator>
    void argsort(
        IndexIterator indices_begin, IndexIterator indices_end,
        ValueIterator data
    ) {
        using Index = typename std::iterator_traits<IndexIterator>::value_type;
        fill_range(indices_begin, indices_end);
        std::sort(indices_begin, indices_end,
        [=](const Index &a, const Index &b) {
            return *(data + a) < *(data + b);
        });
    }
    
    /*!
    \brief Get elements from a sequence according to a sequence of indices.
    
    Let the sequence [`indices_begin`, `indices_end`) be `{i0, i1, ..., iN}`.
    Then the sequence starting at `sorted_data` is filled with
    `{*(data + i0), *(data + i1), ..., *(data + iN)}`.
    */
    template<
        typename IndexIterator,
        typename InputValueIterator,
        typename OutputValueIterator
    >
    void apply(
        IndexIterator indices_begin, IndexIterator indices_end,
        InputValueIterator data,
        OutputValueIterator sorted_data
    ) {
        for (; indices_begin != indices_end; ++indices_begin) {
            *sorted_data = *(data + (*indices_begin));
            ++sorted_data;
        }
    }
}

#endif /* end of include guard: ARGSORT_HPP_FF5D00F6 */
