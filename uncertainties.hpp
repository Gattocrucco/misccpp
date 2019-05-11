#ifndef UNCERTAINTIES_HPP_07A47EC2
#define UNCERTAINTIES_HPP_07A47EC2

#include <map>
#include <numeric>
#include <cmath>
#include <string>
#include <type_traits>
#include <functional>
#include <limits>
#include <cassert>
#include <stdexcept>
#include <ostream>

namespace uncertainties {
    using Id = int;
    
    extern Id next_id;
    
    enum class Order {
        row_major,
        col_major
    };
    
    template<typename Real>
    int _ndigits(const Real &x, const float n) {
        const float log10x = static_cast<float>(std::log10(std::abs(x)));
        const int n_int = static_cast<int>(std::floor(n));
        const float n_frac = n - n_int;
        const float log10x_frac = log10x - std::floor(log10x);
        return n_int + (log10x_frac < n_frac ? 1 : 0);
    }
    
    template<typename Real>
    int _exponent(const Real &x) {
        return static_cast<int>(std::floor(std::log10(std::abs(x))));
    }
    
    template<typename Real>
    std::string _mantissa(const Real &x, const int n, int *const e) {
        const long long m = static_cast<long long>(std::round(x * std::pow(Real(10), n - 1 - *e)));
        std::string s = std::to_string(std::abs(m));
        assert(s.size() == n or s.size() == n + 1 or (m == 0 and n < 0));
        if (n >= 1 and s.size() == n + 1) {
            *e += 1;
            s.pop_back();
        }
        return s;
    }
    
    void _insert_dot(std::string *s, int n, int e) {
        e += s->size() - n;
        n = s->size();
        if (e >= n - 1) {
            // no dot at end of mantissa
        } else if (e >= 0) {
            s->insert(1 + e, 1, '.');
        } else if (e <= -1) {
            s->insert(0, -e, '0');
            s->insert(1, 1, '.');
        }
    }
    
    template<typename Real>
    class UReal {
    private:
        using Type = UReal<Real>;
        
        Real mu;
        std::map<Id, Real> sigma;
        
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
        
    public:
        using real_type = Real;
        
        UReal(const Real n, const Real s):
        mu {std::move(n)}, sigma {{next_id, std::move(s)}} {
            if (s < 0) {
                throw std::invalid_argument("uncertainties::UReal::UReal: s < 0");
            }
            ++next_id;
        }
        
        UReal(const Real &n): UReal(n, 0) {
            ;
        }
        
        const Real &n() const noexcept {
            return this->mu;
        }
        
        Real s() const {
            return std::sqrt(this->s2());
        }
        
        std::string format(const float errdig=1.5f, const std::string &sep=" ± ") const {
            if (errdig <= 1.0f) {
                throw std::invalid_argument("uncertainties::UReal::format: errdig <= 1.0");
            }
            Real s = this->s();
            if (s == 0) {
                return std::to_string(this->mu) + sep + "0";
            }
            const int sndig = _ndigits(s, errdig);
            int sexp = _exponent(s);
            int muexp = this->mu != 0 ? _exponent(this->mu) : sexp - sndig - 1;
            std::string smant = _mantissa(s, sndig, &sexp);
            const int mundig = sndig + muexp - sexp;
            std::string mumant = _mantissa(mu, mundig, &muexp);
            bool use_exp;
            int base_exp;
            if (mundig >= sndig) {
                use_exp = muexp >= mundig or muexp < -1;
                base_exp = muexp;
            } else {
                use_exp = sexp >= sndig or sexp < -1;
                base_exp = sexp;
            }
            if (use_exp) {
                _insert_dot(&mumant, mundig, muexp - base_exp);
                _insert_dot(&smant, sndig, sexp - base_exp);
                return "(" + mumant + sep + smant + ")e" + std::to_string(base_exp);
            } else {
                _insert_dot(&mumant, mundig, muexp);
                _insert_dot(&smant, sndig, sexp);
                return mumant + sep + smant;
            }
        }
        
        operator std::string() {
            return std::to_string(n()) + "+/-" + std::to_string(s());
        }
        
        friend Type copy_unc(const Real n, const Type &x) {
            Type y;
            y.mu = std::move(n);
            y.sigma = x.sigma;
            return y;
        }
        
        template<typename OtherReal>
        friend class UReal;
        
        template<typename OtherReal>
        operator UReal<OtherReal>() {
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
        
        friend Type unary(const Type &x, const Real mu, const Real &dx) {
            Type y;
            y.mu = std::move(mu);
            y.sigma = x.sigma;
            for (auto &it : y.sigma) {
                it.second *= dx;
            }
            return y;
        }
    
        friend Type binary(const Type &x, const Type &y,
                           const Real mu,
                           const Real &dx, const Real &dy) {
            Type z;
            z.mu = std::move(mu);
            for (const auto &it : x.sigma) {
                z.sigma[it.first] = dx * it.second;
            }
            for (const auto &it : y.sigma) {
                z.sigma[it.first] += dy * it.second;
            }
            return z;
        }
        
        template<typename XIt, typename DxIt>
        friend Type nary(XIt xbegin, XIt xend, const Real mu, DxIt dxbegin) {
            Type z;
            z.mu = std::move(mu);
            for (; xbegin != xend; ++xbegin, ++dxbegin) {
                const Type &x = *xbegin;
                const Real &dx = *dxbegin;
                for (const auto &it: x.sigma) {
                    z.sigma[it.first] += dx * it.second;
                }
            }
            return z;
        }
                
        const Type &binary_assign(const Type &x, const Real mu,
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
            this->mu = std::move(mu); // keep this last in case &x == this
            return *this;
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
        friend Type abs(const Type &x) {
            return unary(x, std::abs(x.mu), x.mu >= 0 ? 1 : -1);
        }
        friend Type fmod(const Type &x, const Type &y) {
            return binary(x, y, std::fmod(x.mu, y.mu), 1, -std::trunc(x.mu / y.mu));
        }
        friend Type remainder(const Type &x, const Type &y) {
            return binary(x, y, std::remainder(x.mu, y.mu), 1, -std::round(x.mu / y.mu));
        }
        friend Type fmax(const Type &x, const Type &y) {
            const Real max = std::fmax(x.mu, y.mu);
            const bool c = max == x;
            return binary(x, y, max, c ? 1 : 0, c ? 0 : 1);
        }
        friend Type fmin(const Type &x, const Type &y) {
            const Real min = std::fmin(x.mu, y.mu);
            const bool c = min == x;
            return binary(x, y, min, c ? 1 : 0, c ? 0 : 1);
        }
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
        friend Type erf(const Type &x) {
            static const Real erf_coeff = Real(2) / std::sqrt(Real(3.141592653589793238462643383279502884L));
            return unary(x, std::erf(x.mu), erf_coeff * std::exp(-x.mu * x.mu));
        }
        friend Type erfc(const Type &x) {
            static const Real erf_coeff = Real(2) / std::sqrt(Real(3.141592653589793238462643383279502884L));
            return unary(x, std::erfc(x.mu), -erf_coeff * std::exp(-x.mu * x.mu));
        }
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
    
    template<typename Real>
    std::function<UReal<Real>(const UReal<Real> &)>
    uunary(const std::function<Real(Real)> &f,
           const std::function<Real(Real)> &df) {
        return [f, df](const UReal<Real> &x) {
            return unary(x, f(x.n()), df(x.n()));
        };
    }
    
    template<typename Real>
    constexpr Real default_step() {
        return (1 << (std::numeric_limits<Real>::digits / 2))
               * std::numeric_limits<Real>::epsilon();
    };
    
    template<typename Real>
    std::function<UReal<Real>(const UReal<Real> &)>
    uunary(const std::function<Real(Real)> &f,
           const Real &step=default_step<Real>()) {
        return [f, step](const UReal<Real> &x) {
            const Real &mu = x.n();
            const Real fmu = f(mu);
            const Real dx = (f(mu + step) - fmu) / step;
            return unary(x, fmu, dx);
        };
    }
    
    template<typename Real, typename CharT>
    std::basic_ostream<CharT> &operator<<(std::basic_ostream<CharT> &stream, const UReal<Real> &x) {
        return stream << x.format();
    }
    
    using udouble = UReal<double>;
    using ufloat = UReal<float>;
    
#ifdef UNCERTAINTIES_IMPL
    Id next_id {};
#endif
}

#endif /* end of include guard: UNCERTAINTIES_HPP_07A47EC2 */
