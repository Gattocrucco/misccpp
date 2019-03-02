#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "LinearFit.hpp"

int main() {
    // Generate data.
    std::vector<double> x {1, 2, 3, 4, 5};
    std::vector<double> y(x.size());
    std::transform(x.begin(), x.end(), y.begin(), [](double x) {
        return std::sin(x) + double(std::rand()) / RAND_MAX;
    });
    
    // Fit with LinearFit p[0] + p[1] * sin(x)
    std::vector<double> ones(x.size(), 1.0);
    std::vector<double> sin(x.size());
    std::transform(x.begin(), x.end(), sin.begin(), [](double x) {
        return std::sin(x);
    });
    LinearFit<double> fit(y.begin(), y.end(), ones.begin(), sin.begin());
    std::cout << "LinearFit: \n"
        << fit[0] << "\n"
        << fit[1] << "\nchi2: "
        << fit.chi2() << "\n";
    
    // Again using different constructor
    std::vector<std::vector<double>> fx {ones, sin};
    LinearFit<double> fit2(y, fx);
    std::cout << "LinearFit (2): \n"
        << fit2[0] << "\n"
        << fit2[1] << "\nchi2: "
        << fit2.chi2() << "\n";
    
    // Again
    LinearFit<double> fit3(y, x, {
        [](double x) { return 1; },
        [](double x) { return std::sin(x); }
    });
    std::cout << "LinearFit (3): \n"
        << fit3[0] << "\n"
        << fit3[1] << "\nchi2: "
        << fit3.chi2() << "\n";
}
