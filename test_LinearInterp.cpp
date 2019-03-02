#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "LinearInterp.hpp"

using namespace std;

template<typename Operation>
bool throws(Operation operation) {
    bool throwed;
    try {
        operation();
        throwed = false;
    } catch(exception &e) {
        throwed = true;
    } catch(...) {
        cerr << "throws(): cannot catch exception";
        throw;
    }
    return throwed;
}

int main() {
    const int len {20};
    vector<double> x(len), y(len);
    for (int i = 0; i < len; ++i) {
        x[i] = double(rand()) / RAND_MAX;
        y[i] = double(rand()) / RAND_MAX;
    }
    
    LinearInterp<double, double> interp;
    interp.set_data(x.begin(), x.end(), y.begin());
    
    const int plot_len {1000};
    vector<double> plot_x(plot_len), plot_y(plot_len);
    const auto minmax = minmax_element(x.begin(), x.end());
    const double min {*minmax.first - 0.5}, max {*minmax.second + 0.5};
    for (int i = 0; i < plot_len; ++i) {
        plot_x[i] = min + (max - min) / (plot_len - 1) * i;
        plot_y[i] = interp(
            plot_x[i], decltype(interp)::OutOfRange::Extrapolate
        );
    }
    
    assert(throws([&]() { interp(min - 1); }));
}
