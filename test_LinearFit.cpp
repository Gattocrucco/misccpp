#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cassert>

#include "LinearFit.hpp"

bool close(
    const double x1, const double x2,
    const double atol=1e-8, const double rtol=1e-8
) {
    return std::abs(x1 - x2) <= atol + rtol * (std::abs(x1) + std::abs(x2)) / 2;
}

template<typename Iterator1, typename Iterator2>
bool all_close(
    Iterator1 begin1, Iterator1 end1,
    Iterator2 begin2,
    const double atol=1e-8, const double rtol=1e-8
) {
    for (; begin1 != end1; ++begin1) {
        const double x1 = *begin1, x2 = *begin2;
        if (not close(x1, x2, atol, rtol)) {
            return false;
        }
        ++begin2;
    }
    return true;
}

int main() {
    constexpr double pi = 3.141592653589793;
    
    const std::vector<double> x {0, 1, 2, 3, 4};
    const std::vector<double> y {-1, 0, -1, -2, -1};
    const std::vector<double> expected_parameters {-1, 1};

    LinearFit<double> fit(y, x, { 
        [](double x) { return 1; },
        [](double x) { return std::sin(x * pi/2); }
    });
    
    const std::vector<double> ones(x.size(), 1.0);
    std::vector<double> sinpi2(x.size());
    std::transform(x.begin(), x.end(), sinpi2.begin(), [](double x) {
        return std::sin(x * pi/2);
    });
    
    LinearFit<double> fit2(y, std::vector<std::vector<double>>{ones, sinpi2});
    LinearFit<double> fit3(y.begin(), y.end(), ones.begin(), sinpi2.begin());
    
    for (auto &f : {fit, fit2, fit3}) {
        assert(all_close(f.begin(), f.end(), expected_parameters.begin()));
        assert(close(f.chi2(), 0));
    }
}

// int main() {
//     // Generate data.
//     std::vector<double> x {1, 2, 3, 4, 5};
//     std::vector<double> y(x.size());
//     std::transform(x.begin(), x.end(), y.begin(), [](double x) {
//         return std::sin(x) + double(std::rand()) / RAND_MAX;
//     });
//
//     // Fit with LinearFit p[0] + p[1] * sin(x)
//     std::vector<double> ones(x.size(), 1.0);
//     std::vector<double> sin(x.size());
//     std::transform(x.begin(), x.end(), sin.begin(), [](double x) {
//         return std::sin(x);
//     });
//     LinearFit<double> fit(y.begin(), y.end(), ones.begin(), sin.begin());
//     std::cout << "LinearFit: \n"
//         << fit[0] << "\n"
//         << fit[1] << "\nchi2: "
//         << fit.chi2() << "\n";
//
//     // Again using different constructor
//     std::vector<std::vector<double>> fx {ones, sin};
//     LinearFit<double> fit2(y, fx);
//     std::cout << "LinearFit (2): \n"
//         << fit2[0] << "\n"
//         << fit2[1] << "\nchi2: "
//         << fit2.chi2() << "\n";
//
//     // Again
//     LinearFit<double> fit3(y, x, {
//         [](double x) { return 1; },
//         [](double x) { return std::sin(x); }
//     });
//     std::cout << "LinearFit (3): \n"
//         << fit3[0] << "\n"
//         << fit3[1] << "\nchi2: "
//         << fit3.chi2() << "\n";
// }
