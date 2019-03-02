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
        for (auto fs_it = fs.begin(); fs_it != fs.end(); ++fs_it) {
            size_t row = 0;
            for (auto x_it = std::begin(x); x_it != std::end(x) and row < Ylen; ++x_it) {
                H(row, col) = (*fs_it)(*x_it);
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
    template<
        typename YIterator, typename XIterator, typename... XIterators,
        typename std::enable_if<HasOperatorStar<XIterator>::value, int>::type=0
    >
    explicit LinearFit(
        YIterator ybegin, YIterator yend, XIterator xbegin, XIterators... xs
    ) {
        // explicit first xbegin disambiguates second constructor
        // enable_if disambiguates third constructor
        Y_type Y = Y_builder(ybegin, yend);
        H_type H = H_builder_iterators(Y.size(), xbegin, xs...);
        fit(H, Y);
    }
    
    template<typename YContainer, typename XContainerContainer>
    explicit LinearFit(const YContainer &y, const XContainerContainer &fx) {
        Y_type Y = Y_builder(y);
        H_type H = H_builder_container(Y.size(), fx);
        fit(H, Y);
    }
    
    template<typename YContainer, typename XContainer>
    explicit LinearFit(
        const YContainer &y, const XContainer &x,
        const std::initializer_list<std::function<Scalar(Scalar)>> &fs
    ) {
        Y_type Y = Y_builder(y);
        H_type H = H_builder_functions(Y.size(), x, fs);
        fit(H, Y);
    }
    
    size_t size() const {
        return parameters.size();
    }
    
    Scalar operator[](size_t index) const {
        return parameters(index);
    }
    
    Scalar chi2() const {
        return chi_square;
    }
    
    class Iterator {
    private:
        size_t i;
        const LinearFit<Scalar> * const fit;
    public:
        explicit Iterator(size_t idx, const LinearFit<Scalar> *fit):
            i {idx}, fit {fit}
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
                    "LinearFit<>::Iterator::operator==(): comparing iterators"
                    "from different fit objects"
                );
            }
            return i1.i == i2.i;
        }
        friend bool operator!=(const Iterator &i1, const Iterator &i2) {
            return not (i1 == i2);
        }
    };
    
    Iterator begin() const {
        return Iterator(0, this);
    }
    
    Iterator end() const {
        return Iterator(size(), this);
    }
};

#endif /* end of include guard: LINEARFIT_HPP_BB6921EE */
