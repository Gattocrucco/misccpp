#include <vector>
#include <iostream>
#include <tuple>
#include <type_traits>

#include "tuple_iterator.hpp"

using namespace tuple_iterator;

int main() {
    const std::vector<std::tuple<int, double, unsigned>> v {
        {1, 2.3, 99},
        {3, 4.4, 99},
        {5, 6.0, 99},
        {7, 8.0, 99}
    };
    auto begin1 = tuple_iterator<0>(v.begin());
    auto end1 = tuple_iterator<0>(v.end());
    auto begin2 = tuple_iterator<1>(v.begin());
    auto begin3 = tuple_iterator<2>(v.begin());
    
    std::tuple<int, float> a;
    static_assert(std::is_same<decltype(std::get<0>(a)), int &>::value, "");
    
    while (begin1 != end1) {
        std::cout << *begin1 << ", " << *begin2 << ", " << *begin3 << "\n";
        ++begin1;
        ++begin2;
        ++begin3;
    }
    
    std::vector<std::pair<int, float>> u {
        {3, 7},
        {18, 21}
    };
    auto beginf = tuple_iterator<0>(u.begin());
    auto endf = tuple_iterator<0>(u.end());
    for (auto i = beginf; i != endf; ++i) {
        std::cout << *i << "\n";
        *i = 666;
    }
    for (auto i = beginf; i != endf; ++i) {
        std::cout << *i << "\n";
    }
    
    return 0;
}
