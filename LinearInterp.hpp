#ifndef LINEARINTERP_HPP_F3042564
#define LINEARINTERP_HPP_F3042564

#include <algorithm>
#include <vector>
#include <stdexcept>

#include "argsort.hpp"

/*! \file
\brief Defines class template `LinearInterp` to interpolate linearly a 1D
function.
*/

/*!
\brief Linear 1D interpolator.

The class template parameters are the types of the x and y coordinates. In the
intepolation multiplication is done before division so you can use integers for
both types.

A default constructed interpolator has no data and will throw an exception if it
is used. To set the data, use the `set_data` member function. Then call the
interpolator object on a x coordinate to get the corresponding y.

By default, when computing the interpolated function at a point outside the
range of the data, an exception is thrown. You can change this behavior and
extrapolate from the rightmost or leftmost line by explicitly specifying a flag
each time you call the interpolator.

Example:
~~~{.cpp}
#include <vector>
#include "LinearInterp.hpp"
int main() {
    std::vector<int> x {10, 20, 30, 40, 50};
    std::vector<int> y {0, 1, 2, 1, 0};

    LinearInterp<int, double> interp;
    interp.set_data(x.begin(), x.end(), y.begin(), true);
    // We passed as y a vector of `int`s, but since the second template
    // parameter is `double` the ys will be converted and the output type will
    // be `double`. We use `assume_sorted=true` since we know the `x` vector is
    // already sorted (default is `false`).
    
    double y1 = interp(15); // y1 = 0.5
    // Note: since the x type is `int`, we can only interpolate integer values!
    // If we passed, say, 15.9, it would be silently casted to `int` and
    // yield 15.
    
    double y2 = interp(5, decltype(interp)::OutOfRange::Extrapolate);
    // Since 5 < 10 we are outside the range of the data and we have to pass
    // that wordy flag to extrapolate the value. Now y2 = -0.5.
    
    double y3 = interp(5); // This will throw an exception!
}
~~~
*/
template<typename XValue, typename YValue>
class LinearInterp {
private:
    std::vector<XValue> x;
    std::vector<YValue> y;

public:
    using x_type = XValue;
    using y_type = YValue;
    
    /*!
    \brief Call this function to set the points to be interpolated.
    
    The sequences [`xbegin`, `xend`) and `ybegin` (assumed to have at least the
    same length as the x sequence) are the x and y coordinates of the set of
    points to be interpolated.
    
    You can set `assume_sorted` to `true` if you know that the x sequence is
    sorted.
    
    Note: the x and y values will be casted to the types XValue and YValue, but
    the sorting of the x sequence will take place *before* casting to XValue.
    This should make no difference unless you are using weird user-defined
    types.
    */
    template<typename XIterator, typename YIterator>
    void set_data(
        XIterator xbegin, XIterator xend,
        YIterator ybegin,
        const bool assume_sorted=false
    ) {
        const size_t length = std::distance(xbegin, xend);
        if (length < 2) {
            throw std::runtime_error(
                "LinearInterp::set_data: at least 2 points needed"
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
    
    /*!
    \brief Call the interpolator object to get the y interpolated value
    corresponding to an x.
    
    If `x_point` is outside the range of the points given to `set_data`, the
    behaviour depends on `pl`. If `pl` is `OutOfRange::Throw`,
    `std::runtime_error` is thrown, if it is `OutOfRange::Extrapolate`, the y
    value is obtained using the closer x points.
    */
    YValue operator()(
        const XValue &x_point,
        OutOfRange pl=OutOfRange::Throw
    ) const {
        if (x.size() < 2) {
            throw std::runtime_error(
                "LinearInterp::operator(): at least two points needed,"
                "have you called `set_data` on the interpolator?"
            );
        }
        if (pl == OutOfRange::Throw) {
            if (x_point < x.front() or x.back() < x_point) {
                throw std::runtime_error(
                    "LinearInterp::operator(): x_point outside range"
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
        const XValue &right_x {x[i]};
        const YValue &right_y {y[i]};
        const XValue &left_x {x[i - 1]};
        const YValue &left_y {y[i - 1]};
        if (right_x == left_x) {
            throw std::runtime_error(
                "LinearInterp::operator(): two equal consecutive x in data"
            );
        }
        return left_y +
            (x_point - left_x) * (right_y - left_y) / (right_x - left_x);
    }
};

#endif /* end of include guard: LINEARINTERP_HPP_F3042564 */
