#ifndef LINEARFIT_HPP_BB6921EE
#define LINEARFIT_HPP_BB6921EE

#include <cstddef>
#include <iterator>
#include <cassert>
#include <stdexcept>
#include <string>
#include <initializer_list>
#include <functional>
#include <type_traits>
#include <cmath>

#include <Eigen/Dense>

/*! \file
\brief Defines class `LinearFit` to make linear least-squares fits.

The template library [Eigen](https://eigen.tuxfamily.org) is required. It is a
header only library so it is quick to install.
*/

/*!
\brief Least squares fits of functions linear in the parameters.

This is a nice interface to linear least squares fitting. See
<http://eigen.tuxfamily.org/dox/group__LeastSquares.html> for how to implement
linear least squares fits.

This is "vanilla" least squares, i.e. errors on y and/or x are not supported.

Basic example:
~~~{.cpp}
#include <vector>
#include <cmath>
#include "LinearFit.hpp"
int main() {
    constexpr double pi = 3.141592653589793;
    
    std::vector<double> x {0, 1, 2, 3, 4};
    std::vector<double> y {-1, 0, -1, -2, -1};

    LinearFit<double> fit(y, x, {                   // y =
        [](double x) { return 1; },                 // = A * 1 + 
        [](double x) { return std::sin(x * pi/2); } // + B * sin(x * pi/2)
    }); // NOTE: first y then x
    
    double A = fit[0]; // -1
    double B = fit[1]; // 1
    double chi2 = fit.chi2(); // 0
}
~~~ 
*/
template<typename Scalar>
class LinearFit {
private:
    using size_t = std::size_t;
    using H_type = Eigen::Matrix<
        Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor
    >;
    using Y_type = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    
    void H_builder_iterators_helper(const size_t col, const H_type &H) {
        assert(H.cols() == col);
    }
    
    template<typename XIterator, typename... XIterators>
    void H_builder_iterators_helper(
        const size_t col, H_type &H, XIterator xbegin, XIterators... xs
    ) {
        for (size_t i = 0; i < H.rows(); ++i) {
            H(i, col) = *xbegin;
            ++xbegin;
        }
        H_builder_iterators_helper(col + 1, H, xs...);
    }
    
    template<typename... XIterators>
    H_type H_builder_iterators(const size_t Ylen, XIterators... xs) {
        constexpr size_t Hcols = sizeof...(XIterators);
        H_type H(Ylen, Hcols);
        H_builder_iterators_helper(0, H, xs...);
        return H;
    }
    
    template<typename XContainerContainer>
    H_type H_builder_container(const size_t Ylen, const XContainerContainer &fx) {
        auto fxbegin = std::begin(fx);
        auto fxend = std::end(fx);
        size_t n_fx = std::distance(fxbegin, fxend);
        if (n_fx == 0) {
            throw std::runtime_error("LinearFit: number of parameters is zero");
        }
        H_type H(Ylen, n_fx);
        for (size_t col = 0; fxbegin != fxend; ++fxbegin) {
            auto xbegin = std::begin(*fxbegin);
            auto xend = std::end(*fxbegin);
            size_t row;
            for (row = 0; xbegin != xend and row < Ylen; ++xbegin) {
                H(row, col) = *xbegin;
                ++row;
            }
            if (row < Ylen) {
                throw std::runtime_error(
                    "LinearFit: size of y is " + std::to_string(Ylen) +
                    " but f(x) for parameter " + std::to_string(col) +
                    " has size " + std::to_string(row)
                );
            }
            ++col;
        }
        return H;
    }
    
    template<typename XContainer>
    H_type H_builder_functions(
        const size_t Ylen, const XContainer &x,
        const std::initializer_list<std::function<Scalar(Scalar)>> &fs
    ) {
        H_type H(Ylen, fs.size());
        size_t col = 0;
        for (const auto &f : fs) {
            size_t row = 0;
            for (auto x_it = std::begin(x); x_it != std::end(x) and row < Ylen; ++x_it) {
                H(row, col) = f(*x_it);
                ++row;
            }
            if (row < Ylen) {
                throw std::runtime_error(
                    "LinearFit: size of y is " + std::to_string(Ylen) +
                    " but f(x) for parameter " + std::to_string(col) +
                    " has size " + std::to_string(row)
                );
            }
            ++col;
        }
        return H;
    }

    template<typename YIterator>
    Y_type Y_builder(YIterator begin, YIterator end) {
        const size_t len = std::distance(begin, end);
        if (len == 0) {
            throw std::runtime_error("LinearFit: Y has zero elements");
        }
        Y_type Y(len);
        for (size_t i = 0; i < len; ++i) {
            Y(i) = *begin;
            ++begin;
        }
        return Y;
    }
    
    template<typename YContainer>
    Y_type Y_builder(const YContainer &y) {
        return Y_builder(std::begin(y), std::end(y));
    }
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> parameters;
    Scalar chi_square;
    
    void fit_normal_equations(const H_type &H, const Y_type &Y) {
        // parameters = (H^T * H)^(-1) * H^T * Y
        // this is the quickest and easier to understand method,
        // but use `fit_svd` instead
        auto HT = H.transpose();
        auto HTH = HT * H;
        auto HTY = HT * Y;
        parameters = HTH.colPivHouseholderQr().solve(HTY);
        chi_square = (Y - H * parameters).squaredNorm();
    }
    
    void fit_svd(const H_type &H, const Y_type &Y) {
        parameters = H.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Y);
        chi_square = (Y - H * parameters).squaredNorm();
    }
    
    void fit(const H_type &H, const Y_type &Y) {
        fit_svd(H, Y);
        if (std::isnan(chi_square)) {
            throw std::runtime_error("LinearFit: NaN chi2");
        }
    }
    
    template<typename T>
    class HasOperatorStar {
    private:
        using Yes = char[1];
        using No = char[2];
        template<typename C> static Yes &test(decltype(&C::operator*));
        template<typename C> static No &test(...);
    public:
        constexpr static bool value = sizeof(test<T>(0)) == sizeof(Yes);
    };
    
public:
    /*!
    \brief Fit y = p0 * x0 + p1 * x1 + ...
    
    The y, x0, x1, ... sequences are passed by iterators, so this is the most
    general constructor. Note that x0, x1, etc. may not be what you call "x" in
    your problem, tipically you will call them f0(x), f1(x), etc.
    
    The constructor will throw `std::runtime_error` if the chisquare ends up
    NaN or if the y sequence is empty.
    */
    template<
        typename YIterator, typename XIterator, typename... XIterators,
        typename std::enable_if<HasOperatorStar<XIterator>::value, int>::type=0
    >
    explicit LinearFit(
        YIterator ybegin, YIterator yend, XIterator xbegin, XIterators... xbegins
    ) {
        // explicit first xbegin disambiguates second constructor
        // enable_if disambiguates third constructor
        Y_type Y = Y_builder(ybegin, yend);
        H_type H = H_builder_iterators(Y.size(), xbegin, xbegins...);
        fit(H, Y);
    }
    
    /*!
    \brief Fit y = p0 * fx[0] + p1 * fx[1] + ...
    
    The argument `y` must be a container i.e. something that `std::begin` and
    `std::end` understand. `fx` must be a container of containers.
    
    The constructor will throw `std::runtime_error` if the chisquare ends up
    NaN, if `y` is empty, if `fx` is empty or if any of the containers inside
    `fx` is shorter than `y`.
    */
    template<typename YContainer, typename XContainerContainer>
    explicit LinearFit(const YContainer &y, const XContainerContainer &fx) {
        Y_type Y = Y_builder(y);
        H_type H = H_builder_container(Y.size(), fx);
        fit(H, Y);
    }
    
    /*!
    \brief Fit y = p0 * fs[0](x) + p1 * fs[1](x) + ...
    
    The arguments `y` and `x` must be containers i.e. something that
    `std::begin` and `std::end` understand.
    
    `fs` is a `std::intializer_list` of functions. Typically you will write it
    directly in the function call as a braced list of lambdas:
    ~~~{.cpp}
    LinearFit<Scalar>(y, x, { [](Scalar x) { return x; }, ...});
    ~~~
    
    The constructor will throw `std::runtime_error` if the chisquare ends up
    NaN, if `y` is empty or if `fs` is empty.
    */
    template<typename YContainer, typename XContainer>
    explicit LinearFit(
        const YContainer &y, const XContainer &x,
        const std::initializer_list<std::function<Scalar(Scalar)>> &fs
    ) {
        Y_type Y = Y_builder(y);
        H_type H = H_builder_functions(Y.size(), x, fs);
        fit(H, Y);
    }
    
    /*!
    \brief Returns the number of parameters.
    */
    size_t size() const {
        return parameters.size();
    }
    
    /*!
    \brief Returns the `index`th fitted parameter.
    
    The indices of the parameters are 0-based and are implicitly defined by the
    order of the functions when constructing the fit object.
    */
    Scalar operator[](size_t index) const {
        return parameters(index);
    }
    
    /*!
    \brief Returns the chisquared.
    */
    Scalar chi2() const {
        return chi_square;
    }
    
    class Iterator {
    private:
        size_t i;
        const LinearFit<Scalar> * const fit;
    public:
        explicit Iterator(size_t idx, const LinearFit<Scalar> *fit_obj):
            i {idx}, fit {fit_obj}
        { ; }
        void operator++() {
            ++i;
        }
        Scalar operator*() const {
            return (*fit)[i];
        }
        friend bool operator==(const Iterator &i1, const Iterator &i2) {
            if (i1.fit != i2.fit) {
                throw std::runtime_error(
                    "LinearFit::Iterator::operator==: comparing iterators "
                    "from different fit objects"
                );
            }
            return i1.i == i2.i;
        }
        friend bool operator!=(const Iterator &i1, const Iterator &i2) {
            return not (i1 == i2);
        }
    };
    
    /*!
    \brief Returns an iterator to the parameters pointing at the first one.
    */
    Iterator begin() const {
        return Iterator(0, this);
    }
    
    /*!
    \brief Returns an iterator to the parameters pointing at one past the last
    one.
    */
    Iterator end() const {
        return Iterator(size(), this);
    }
};

#endif /* end of include guard: LINEARFIT_HPP_BB6921EE */
