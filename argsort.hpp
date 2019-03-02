#ifndef ARGSORT_HPP_FF5D00F6
#define ARGSORT_HPP_FF5D00F6

#include <algorithm>
#include <iterator>

namespace argsort {
    template<typename Iterator>
    void fill_range(Iterator begin, Iterator end) {
        using Type = typename std::iterator_traits<Iterator>::value_type;
        for (Type i {}; begin != end; ++begin) {
            *begin = i;
            ++i;
        }
    }
    
    template<typename IndexIterator, typename ValueIterator>
    void argsort(
        IndexIterator indices_begin, IndexIterator indices_end,
        ValueIterator data
    ) {
        using Index = typename std::iterator_traits<IndexIterator>::value_type;
        fill_range(indices_begin, indices_end);
        std::sort(indices_begin, indices_end, [=](Index a, Index b) {
            return *(data + a) < *(data + b);
        });
    }
    
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
        for (IndexIterator idx = indices_begin; idx != indices_end; ++idx) {
            *sorted_data = *(data + (*idx));
            ++sorted_data;
        }
    }
}

#endif /* end of include guard: ARGSORT_HPP_FF5D00F6 */
