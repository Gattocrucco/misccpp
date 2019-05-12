#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>

#include <uncertainties/stat.hpp>
#include <uncertainties/io.hpp>
#include <uncertainties/math.hpp>
#include <uncertainties/ureal.hpp>

namespace unc = uncertainties;
using unc::udouble;
using unc::ufloat;

template<typename Vector>
void print_matrix(Vector m) {
    const int n = static_cast<int>(std::round(std::sqrt(m.size())));
    char s[1024];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::snprintf(s, sizeof(s), "%10.6g, ", double(m[n * i + j]));
            std::cout << s;
        }
        std::cout << "\n";
    }
}

double sin2(double x) {
    return std::sin(x) * std::sin(x);
}

double dsin2(double x) {
    return 2 * std::sin(x) * std::cos(x);
}

int main() {
    udouble x(1, 1);
    udouble y {1, 1};
    std::cout << "x - x = " << (x - x).format() << "\n";
    std::cout << "x + x = " << (x + x).format() << "\n";
    std::cout << "x - y = " << (x - y).format() << "\n";
    std::cout << "sin(x) = " << sin(x).format() << "\n";
    std::cout << "sizeof(udouble) = " << sizeof(udouble) << "\n";
    std::cout << "ufloat(x) = " << ufloat(x).format() << "\n";
    std::cout << "x - ufloat(x) = " << (x - udouble(ufloat(x))).format() << "\n";
    std::cout << "cov(x, x) = " << cov(x, x) << "\n";
    std::cout << "cov(x, -x) = " << cov(x, -x) << "\n";
    std::cout << "cov(x, y) = " << cov(x, y) << "\n";
    std::cout << "cov(x, x + y) = " << cov(x, x + y) << "\n";
    x += 1;
    std::cout << "x += 1; x = " << x.format() << "\n";
    x += x;
    std::cout << "x += x; x = " << x.format() << "\n";
    std::cout << "x * x = " << (x * x).format() << "\n";
    x *= x;
    std::cout << "x *= x; x = " << x.format() << "\n";
    std::vector<double> corr = unc::corr_matrix<std::vector<double>>(std::vector<udouble>{x, 2 * x, y, x + y});
    std::cout << "corr_matrix({x, 2 * x, y, x + y}) = \n";
    print_matrix(corr);
    std::vector<udouble> v {x, 2 * x, y, x + y};
    std::vector<double> cov(v.size() * v.size());
    unc::cov_matrix(v.begin(), v.end(), cov.begin());
    std::cout << "cov_matrix({x, 2 * x, y, x + y}) = \n";
    print_matrix(cov);
    auto usin2 = unc::uunary<double>(sin2, dsin2);
    auto usin2num = unc::uunary<double>(sin2);
    std::cout << "sin2(y) = " << usin2(y).format(15) << "\n";
    std::cout << "sin2(y) (auto deriv) = " << usin2num(y).format(15) << "\n";
    std::cout << udouble(0, 0.001).format(2) << "\n";
    std::cout << udouble(0, 0.1).format(2) << "\n";
    std::cout << udouble(0, 1).format(2) << "\n";
    std::cout << udouble(0, 10).format(2) << "\n";
    std::cout << udouble(0, 100).format(2) << "\n";
    std::cout << udouble(0, 0.0196).format(2) << "\n";
    std::cout << udouble(0, 0.196).format(2) << "\n";
    std::cout << udouble(0, 1.96).format(2) << "\n";
    std::cout << udouble(0, 19.6).format(2) << "\n";
    std::cout << udouble(0, 196).format(2) << "\n";
    std::cout << udouble(0, 0.0996).format(2) << "\n";
    std::cout << udouble(0, 0.996).format(2) << "\n";
    std::cout << udouble(0, 9.96).format(2) << "\n";
    std::cout << udouble(0, 99.6).format(2) << "\n";
    std::cout << udouble(0, 996).format(2) << "\n";
    std::cout << udouble(0.025, 0.003).format(2) << "\n";
    std::cout << udouble(0.025, 0.0003).format(2) << "\n";
    std::cout << udouble(0.025, 0.00003).format(2) << "\n";
    std::cout << udouble(0.003  , 0.025).format(2) << "\n";
    std::cout << udouble(0.0003 , 0.025).format(2) << "\n";
    std::cout << udouble(0.00003, 0.025).format(2) << "\n";
    std::cout << udouble(0.25, 0.3).format(2) << "\n";
    std::cout << udouble(0.25, 0.03).format(2) << "\n";
    std::cout << udouble(0.25, 0.003).format(2) << "\n";
    std::cout << udouble(0.25, 0.0003).format(2) << "\n";
    std::cout << udouble(-1, 1) << "\n";
    std::cout << unc::format(2345678901) << "\n";
}
