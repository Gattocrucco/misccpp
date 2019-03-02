#ifndef LINEARINTERP_HPP_F3042564
#define LINEARINTERP_HPP_F3042564

#include <algorithm>
#include <vector>
#include <stdexcept>

#include "argsort.hpp"

template<typename XValue, typename YValue, typename ComputeType=double>
class LinearInterp {
private:
    std::vector<XValue> x;
    std::vector<YValue> y;

public:
    using x_type = XValue;
    using y_type = YValue;
    using compute_type = ComputeType;
        
    template<typename XIterator, typename YIterator>
    void set_data(
        XIterator xbegin, XIterator xend,
        YIterator ybegin,
        const bool assume_sorted=false
    ) {
        const size_t length = std::distance(xbegin, xend);
        if (length < 2) {
            throw std::runtime_error(
                "LinearInterp: set_data(): at least 2 points needed"
            );
        }
        x.resize(length);
        y.resize(length);
        if (assume_sorted) {
            std::copy(xbegin, xend, x.begin());
            std::copy(ybegin, ybegin + length, y.begin());
        } else {
            std::vector<size_t> indices(length);
            argsort::argsort(indices.begin(), indices.end(), xbegin);
            argsort::apply(indices.begin(), indices.end(), xbegin, x.begin());
            argsort::apply(indices.begin(), indices.end(), ybegin, y.begin());
        }
    }
    
    enum class OutOfRange {
        Throw, Extrapolate
    };
    
    YValue operator()(XValue x_point, OutOfRange pl=OutOfRange::Throw) const {
        if (pl == OutOfRange::Throw) {
            if (x_point < x.front() or x.back() < x_point) {
                throw std::runtime_error(
                    "LinearInterp: operator(): x_point outside range"
                );
            }
        }
        size_t i = std::distance(
            x.begin(), std::lower_bound(x.begin(), x.end(), x_point)
        );
        // `i` is the lowest such that `x_point <= x[i]`. We will use x[i] and
        // x[i - 1] to interpolate.
        if (i == 0) {
            // If we checked the bounds, it implies `x_point == x[0]`. Otherwise
            // anyway we need the first two points to extrapolate.
            i = 1;
        } else if (i == x.size()) {
            // It implies `x_point > x.back()` and that we have to extrapolate.
            i = x.size() - 1;
        }
        const ComputeType right_x {x[i]};
        const ComputeType right_y {y[i]};
        const ComputeType left_x {x[i - 1]};
        const ComputeType left_y {y[i - 1]};
        if (right_x == left_x) {
            throw std::runtime_error(
                "LinearInterp: operator(): two equal consecutive x in data."
            );
        }
        const ComputeType slope {(right_y - left_y) / (right_x - left_x)};
        const YValue y_point = left_y + (ComputeType{x_point} - left_x) * slope;
        return y_point;
    }
};

#endif /* end of include guard: LINEARINTERP_HPP_F3042564 */
