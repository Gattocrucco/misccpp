#ifndef UNCERTAINTIES_STAT_HPP_C71B9D96
#define UNCERTAINTIES_STAT_HPP_C71B9D96

/*! \file
\brief Defines functions to manipulate covariance matrices.
*/

#include <iterator>

#include "core.hpp"

namespace uncertainties {
    namespace internal {
        template<typename InputIt, typename OutputIt, typename Operation>
        OutputIt outer(InputIt begin, InputIt end, OutputIt matrix,
                       Operation op, Order order=Order::row_major) {
            for (InputIt i = begin; i != end; ++i) {
                for (InputIt j = begin; j != end; ++j) {
                    *matrix = order == Order::row_major ? op(*i, *j) : op(*j, *i);
                    ++matrix;
                }
            }
            return matrix;
        }

        template<typename OutVector, typename InVector, typename Operation>
        OutVector outer(InVector x, Operation op, Order order=Order::row_major) {
            const typename InVector::size_type n = x.size();
            OutVector matrix(n * n);
            internal::outer(std::begin(x), std::end(x), std::begin(matrix), op, order);
            return matrix;
        }
    }
    
    template<typename InputIt, typename OutputIt>
    OutputIt cov_matrix(InputIt begin, InputIt end, OutputIt matrix,
                        Order order=Order::row_major) {
        using Type = typename InputIt::value_type;
        return internal::outer(begin, end, matrix, [](const Type &x, const Type &y) {
            return cov(x, y);
        }, order);
    }
    
    template<typename OutVector, typename InVector>
    OutVector cov_matrix(InVector x, Order order=Order::row_major) {
        using Type = typename InVector::value_type;
        return internal::outer<OutVector>(x, [](const Type &x, const Type &y) {
            return cov(x, y);
        }, order);
    }
    
    template<typename InputIt, typename OutputIt>
    OutputIt corr_matrix(InputIt begin, InputIt end, OutputIt matrix,
                         Order order=Order::row_major) {
        using Type = typename InputIt::value_type;
        return internal::outer(begin, end, matrix, [](const Type &x, const Type &y) {
            return corr(x, y);
        }, order);
    }
    
    template<typename OutVector, typename InVector>
    OutVector corr_matrix(InVector x, Order order=Order::row_major) {
        using Type = typename InVector::value_type;
        return internal::outer<OutVector>(x, [](const Type &x, const Type &y) {
            return corr(x, y);
        }, order);
    }
}

#endif /* end of include guard: UNCERTAINTIES_STAT_HPP_C71B9D96 */
