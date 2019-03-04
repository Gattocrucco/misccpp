#include <iostream>
#include <vector>
#include <stdexcept>

#include "LinearInterp.hpp"

template<typename Operation>
bool throws(Operation operation) {
    bool throwed;
    try {
        operation();
        throwed = false;
    } catch(std::exception &e) {
        throwed = true;
    } catch(...) {
        std::cerr << "throws(): cannot catch exception";
        throw;
    }
    return throwed;
}

int main() {
    std::vector<int> x {10, 20, 30, 40, 50};
    std::vector<int> y {0, 1, 2, 1, 0};
    
    LinearInterp<int, double> interp;
    assert(throws([&]() {
        interp.set_data(x.begin(), x.begin() + 1, y.begin());
    }));
    assert(throws([&]() {
        interp(10);
    }));

    interp.set_data(x.begin(), x.end(), y.begin(), true);
    
    assert(interp(15) == 0.5);
    
    assert(interp(5, decltype(interp)::OutOfRange::Extrapolate) == -0.5);
    
    assert(throws([&]() {
        interp(5);
    }));
}
