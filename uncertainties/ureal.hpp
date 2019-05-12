#ifndef UNCERTAINTIES_UREAL_HPP_07A47EC2
#define UNCERTAINTIES_UREAL_HPP_07A47EC2

/*! \file
\brief Defines class template `UReal` and basic utilities.
*/

#include <map>
#include <cmath>
#include <string>
#include <functional>
#include <limits>
#include <cassert>
#include <stdexcept>

/*!
\brief C++ header library for linear uncertainty propagation.

Basic example:
~~~{.cpp}
#include <iostream>
#include <uncertainties/ureal.hpp>
namespace unc = uncertainties;
int main() {
    unc::udouble x(2, 1), y(2, 1);
    unc::udouble y = x - x;
    unc::udouble z = x - y;
    std::cout << y.format(2) << ", " << z.format(2) << "\n";
}
~~~
*/
namespace uncertainties {
    namespace internal {
        using Id = int;
    
        extern Id next_id;
    
        template<typename Real>
        int ndigits(const Real &x, const float n) {
            const float log10x = static_cast<float>(std::log10(std::abs(x)));
            const int n_int = static_cast<int>(std::floor(n));
            const float n_frac = n - n_int;
            const float log10x_frac = log10x - std::floor(log10x);
            return n_int + (log10x_frac < n_frac ? 1 : 0);
        }
    
        template<typename Real>
        int exponent(const Real &x) {
            return static_cast<int>(std::floor(std::log10(std::abs(x))));
        }
    
        template<typename Real>
        std::string mantissa(const Real &x, const int n, int *const e) {
            const long long m = static_cast<long long>(
                std::round(x * std::pow(Real(10), n - 1 - *e))
            );
            std::string s = std::to_string(std::abs(m));
            assert(s.size() == n or s.size() == n + 1 or (m == 0 and n < 0));
            if (n >= 1 and s.size() == n + 1) {
                *e += 1;
                s.pop_back();
            }
            return s;
        }
    
        void insert_dot(std::string *s, int n, int e);
        std::string format_exp(const int e);
    }
    
    /*!
    \brief Represents a number with associated uncertainty.
    */
    template<typename Real>
    class UReal {
    private:
        using Type = UReal<Real>;
        
        Real mu;
        std::map<internal::Id, Real> sigma;
        
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
        mu {std::move(n)}, sigma {{internal::next_id, std::move(s)}} {
            if (s < 0) {
                throw std::invalid_argument("uncertainties::UReal::UReal: s < 0");
            }
            ++internal::next_id;
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
        
        std::string format(const float errdig=1.5f,
                           const std::string &sep=" ± ") const {
            if (errdig <= 1.0f) {
                throw std::invalid_argument("uncertainties::UReal::format: errdig <= 1.0");
            }
            Real s = this->s();
            if (s == 0) {
                return std::to_string(this->mu) + sep + "0";
            }
            const int sndig = internal::ndigits(s, errdig);
            int sexp = internal::exponent(s);
            int muexp = this->mu != 0 ? internal::exponent(this->mu) : sexp - sndig - 1;
            std::string smant = internal::mantissa(s, sndig, &sexp);
            const int mundig = sndig + muexp - sexp;
            std::string mumant = internal::mantissa(this->mu, mundig, &muexp);
            const std::string musign = this->mu < 0 ? "-" : "";
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
                internal::insert_dot(&mumant, mundig, muexp - base_exp);
                internal::insert_dot(&smant, sndig, sexp - base_exp);
                return "(" + musign + mumant + sep + smant + ")e" + internal::format_exp(base_exp);
            } else {
                internal::insert_dot(&mumant, mundig, muexp);
                internal::insert_dot(&smant, sndig, sexp);
                return musign + mumant + sep + smant;
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
    };
    
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
    
    using udouble = UReal<double>;
    using ufloat = UReal<float>;
}

#endif /* end of include guard: UNCERTAINTIES_UREAL_HPP_07A47EC2 */
