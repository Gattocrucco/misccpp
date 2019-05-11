#ifndef UNCERTAINTIES_HPP_07A47EC2
#define UNCERTAINTIES_HPP_07A47EC2

#include <map>
#include <numeric>
#include <cmath>
#include <string>
#include <type_traits>

namespace uncertainties {
    using ID = int;
    
    extern ID next_id;
    
    class NegativeSigma {};
    
    enum class Order {
        row_major,
        col_major
    };
    
    template<typename Real>
    class UReal {
    private:
        using Type = UReal<Real>;
        
        Real mu;
        std::map<ID, Real> sigma;
        
        UReal() {
            ;
        }
        
        Real s2() const {
            Real s2(0);
            for (const auto &it : this->sigma) {
                const Real &s = it.second;
                s2 += s * s;
            }
            return s2;
        }
        
        friend Type unary(const Type &x, const Real &mu, const Real &dx) {
            Type y;
            y.mu = mu;
            y.sigma = x.sigma;
            for (auto &it : y.sigma) {
                it.second *= dx;
            }
            return y;
        }
    
        friend Type binary(const Type &x, const Type &y,
                           const Real &mu, const Real &dx, const Real &dy) {
            Type z;
            z.mu = mu;
            for (const auto &it : x.sigma) {
                z.sigma[it.first] = dx * it.second;
            }
            for (const auto &it : y.sigma) {
                z.sigma[it.first] += dy * it.second;
            }
            return z;
        }
        
        const Type &binary_assign(const Type &x, const Real &mu,
                                  const Real &dt, const Real &dx) {
            if (&x == this) {
                const Real d = dt + dx;
                for (auto &it : this->sigma) {
                    it.second *= d;
                }
            } else {
                if (dt != 1) {
                    for (auto &it : this->sigma) {
                        it.second *= dt;
                    }
                }
                for (const auto &it : x.sigma) {
                    this->sigma[it.first] += dx * it.second;
                }
            }
            this->mu = mu; // keep this last in case &x == this
            return *this;
        }
        
    public:
        using real_type = Real;
        
        UReal(const Real &n, const Real &s) {
            if (s < 0) {
                throw NegativeSigma();
            }
            this->mu = n;
            this->sigma[next_id] = s;
            ++next_id;
        }
        
        UReal(const Real &n): UReal(n, 0) {
            ;
        }
        
        Real n() const {
            return this->mu;
        }
        
        Real s() const {
            return std::sqrt(this->s2());
        }
        
        std::string format() {
            return std::to_string(n()) + " ± " + std::to_string(s());
        }
        
        friend Type copy_unc(const Real &n, const Type &x) {
            Type y;
            y.mu = n;
            y.sigma = x.sigma;
            return y;
        }
        
        template<typename OtherReal>
        friend class UReal;
        
        template<typename OtherReal>
        operator UReal<OtherReal>() {
            static_assert(!std::is_same<Real, OtherReal>::value, "what??");
            UReal<OtherReal> x;
            x.mu = this->mu;
            for (const auto &it : this->sigma) {
                x.sigma[it.first] = it.second;
            }
            return x;
        }
        
        friend Real cov(const Type &x, const Type &y) {
            Real cov(0);
            const Type *min_size, *max_size;
            if (x.sigma.size() > y.sigma.size()) {
                min_size = &y;
                max_size = &x;
            } else {
                min_size = &x;
                max_size = &y;
            }
            const auto end = max_size->sigma.end();
            for (const auto &it : min_size->sigma) {
                const auto IT = max_size->sigma.find(it.first);
                if (IT != end) {
                    cov += it.second * IT->second;
                }
            }
            return cov;
        }
        
        friend Real var(const Type &x) {
            return x.s2();
        }
        
        friend Real corr(const Type &x, const Type &y) {
            return cov(x, y) / std::sqrt(var(x) * var(y));
        }
        
        friend Type operator+(const Type &x) {
            return unary(x, x.mu, 1);
        }
        friend Type operator-(const Type &x) {
            return unary(x, -x.mu, -1);
        }
        friend Type operator+(const Type &x, const Type &y) {
            return binary(x, y, x.mu + y.mu, 1, 1);
        }
        friend Type operator-(const Type &x, const Type &y) {
            return binary(x, y, x.mu - y.mu, 1, -1);
        }
        friend Type operator*(const Type &x, const Type &y) {
            return binary(x, y, x.mu * y.mu, y.mu, x.mu);
        }
        friend Type operator/(const Type &x, const Type &y) {
            const Real inv_y = Real(1) / y.mu;
            const Real mu = x.mu * inv_y;
            return binary(x, y, mu, inv_y, -mu * inv_y);
        }
        const Type &operator+=(const Type &x) {
            return binary_assign(x, this->mu + x.mu, 1, 1);
        }
        const Type &operator-=(const Type &x) {
            return binary_assign(x, this->mu - x.mu, 1, -1);
        }
        const Type &operator*=(const Type &x) {
            return binary_assign(x, this->mu * x.mu, x.mu, this->mu);
        }
        const Type &operator/=(const Type &x) {
            const Real inv_x = Real(1) / x.mu;
            const Real mu = this->mu * inv_x;
            return binary_assign(x, mu, inv_x, -mu * inv_x);
        }
        // abs?
        // fmod?
        // remainder?
        // fma?
        // fmax?
        // fmin?
        // fdim?
        friend Type exp(const Type &x) {
            return unary(x, std::exp(x.mu), std::exp(x.mu));
        }
        friend Type exp2(const Type &x) {
            return unary(x, std::exp2(x.mu), std::log(Real(2)) * std::exp(x.mu)); 
        }
        friend Type expm1(const Type &x) {
            return unary(x, std::expm1(x.mu), std::exp(x.mu));
        }
        friend Type log(const Type &x) {
            return unary(x, std::log(x.mu), Real(1) / x.mu);
        }
        friend Type log10(const Type &x) {
            return unary(x, std::log10(x.mu), Real(1) / (x.mu * std::log(Real(10))));
        }
        friend Type log2(const Type &x) {
            return unary(x, std::log2(x.mu), Real(1) / (x.mu * std::log(Real(2))));
        }
        friend Type log1p(const Type &x) {
            return unary(x, std::log1p(x.mu), Real(1) / (Real(1) + x.mu));
        }
        friend Type pow(const Type &x, const Type &y) {
            const Real p = std::pow(x.mu, y.mu);
            return binary(x, y, p, p * y.mu / x.mu, p * std::log(x.mu));
        }
        friend Type sqrt(const Type &x) {
            return unary(x, std::sqrt(x.mu), Real(1) / (2 * std::sqrt(x.mu)));
        }
        friend Type cbrt(const Type &x) {
            return unary(x, std::cbrt(x.mu), std::pow(x.mu, -Real(2) / Real(3)) / Real(3));
        }
        friend Type hypot(const Type &x, const Type &y) {
            const Real h = std::hypot(x.mu, y.mu);
            return binary(x, y, h, x.mu / h, y.mu / h);
        }
        friend Type sin(const Type &x) {
            return unary(x, std::sin(x.mu), std::cos(x.mu));
        }
        friend Type cos(const Type &x) {
            return unary(x, std::cos(x.mu), -std::sin(x.mu));
        }
        friend Type tan(const Type &x) {
            const Real t = std::tan(x.mu);
            return unary(x, t, Real(1) + t * t);
        }
        friend Type asin(const Type &x) {
            return unary(x, std::asin(x.mu), Real(1) / std::sqrt(1 - x.mu * x.mu));
        }
        friend Type acos(const Type &x) {
            return unary(x, std::acos(x.mu), -Real(1) / std::sqrt(1 - x.mu * x.mu));
        }
        friend Type atan(const Type &x) {
            return unary(x, std::atan(x.mu), Real(1) / (1 + x.mu * x.mu));
        }
        friend Type atan2(const Type &x, const Type &y) {
            const Real yx = y.mu / x.mu;
            const Real dy = Real(1) / ((Real(1) + yx * yx) * x);
            const Real dx = dy * (-yx);
            return binary(x, y, std::atan2(x.mu, y.mu), dx, dy);
        }
        friend Type sinh(const Type &x) {
            return unary(x, std::sinh(x.mu), std::cosh(x.mu));
        }
        friend Type cosh(const Type &x) {
            return unary(x, std::cosh(x.mu), std::sinh(x.mu));
        }
        friend Type tanh(const Type &x) {
            const Real t = std::tanh(x.mu);
            return unary(x, t, Real(1) - t * t);
        }
        friend Type asinh(const Type &x) {
            return unary(x, std::asinh(x.mu), Real(1) / std::sqrt(x.mu * x.mu + 1));
        }
        friend Type acosh(const Type &x) {
            return unary(x, std::acosh(x.mu), Real(1) / std::sqrt(x.mu * x.mu - 1));
        }
        friend Type atanh(const Type &x) {
            return unary(x, std::atanh(x.mu), Real(1) / (1 - x.mu * x.mu));
        }
        // erf
        // erfc
        // tgamma
        // lgamma
        // copysign?
        friend bool isfinite(const Type &x) {
            return std::isfinite(x.mu) and std::isfinite(x.s2());
        }
        friend bool isnormal(const Type &x) {
            return std::isnormal(x.mu) and std::isnormal(x.s2());
        }
    };
        
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
        outer(std::begin(x), std::end(x), std::begin(matrix), op, order);
        return matrix;
    }
    
    template<typename InputIt, typename OutputIt>
    OutputIt cov_matrix(InputIt begin, InputIt end, OutputIt matrix,
                        Order order=Order::row_major) {
        using Type = typename InputIt::value_type;
        return outer(begin, end, matrix, [](const Type &x, const Type &y) {
            return cov(x, y);
        }, order);
    }
    
    template<typename OutVector, typename InVector>
    OutVector cov_matrix(InVector x, Order order=Order::row_major) {
        using Type = typename InVector::value_type;
        return outer<OutVector>(x, [](const Type &x, const Type &y) {
            return cov(x, y);
        }, order);
    }
    
    template<typename InputIt, typename OutputIt>
    OutputIt corr_matrix(InputIt begin, InputIt end, OutputIt matrix,
                        Order order=Order::row_major) {
        using Type = typename InputIt::value_type;
        return outer(begin, end, matrix, [](const Type &x, const Type &y) {
            return corr(x, y);
        }, order);
    }
    
    template<typename OutVector, typename InVector>
    OutVector corr_matrix(InVector x, Order order=Order::row_major) {
        using Type = typename InVector::value_type;
        return outer<OutVector>(x, [](const Type &x, const Type &y) {
            return corr(x, y);
        }, order);
    }
    
    using udouble = UReal<double>;
    using ufloat = UReal<float>;
    
#ifdef UNCERTAINTIES_IMPL
    ID next_id {};
#endif
}

#endif /* end of include guard: UNCERTAINTIES_HPP_07A47EC2 */
