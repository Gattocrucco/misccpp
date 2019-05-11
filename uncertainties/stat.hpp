#ifndef UNCERTAINTIES_STAT_HPP_C71B9D96
#define UNCERTAINTIES_STAT_HPP_C71B9D96

namespace uncertainties {
    template<typename Real>
    class UReal;

    enum class Order {
        row_major,
        col_major
    };
    
    template<typename InputIt, typename OutputIt, typename Operation>
    OutputIt _outer(InputIt begin, InputIt end, OutputIt matrix,
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
    OutVector _outer(InVector x, Operation op, Order order=Order::row_major) {
        const typename InVector::size_type n = x.size();
        OutVector matrix(n * n);
        _outer(std::begin(x), std::end(x), std::begin(matrix), op, order);
        return matrix;
    }
    
    template<typename InputIt, typename OutputIt>
    OutputIt cov_matrix(InputIt begin, InputIt end, OutputIt matrix,
                        Order order=Order::row_major) {
        using Type = typename InputIt::value_type;
        return _outer(begin, end, matrix, [](const Type &x, const Type &y) {
            return cov(x, y);
        }, order);
    }
    
    template<typename OutVector, typename InVector>
    OutVector cov_matrix(InVector x, Order order=Order::row_major) {
        using Type = typename InVector::value_type;
        return _outer<OutVector>(x, [](const Type &x, const Type &y) {
            return cov(x, y);
        }, order);
    }
    
    template<typename InputIt, typename OutputIt>
    OutputIt corr_matrix(InputIt begin, InputIt end, OutputIt matrix,
                        Order order=Order::row_major) {
        using Type = typename InputIt::value_type;
        return _outer(begin, end, matrix, [](const Type &x, const Type &y) {
            return corr(x, y);
        }, order);
    }
    
    template<typename OutVector, typename InVector>
    OutVector corr_matrix(InVector x, Order order=Order::row_major) {
        using Type = typename InVector::value_type;
        return _outer<OutVector>(x, [](const Type &x, const Type &y) {
            return corr(x, y);
        }, order);
    }
}

#endif /* end of include guard: UNCERTAINTIES_STAT_HPP_C71B9D96 */
