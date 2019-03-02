#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cassert>

#include "argsort.hpp"

int main() {
    // Fill `array` with random numbers.
    std::vector<int> array;
    for (int i = 0; i < 1000; ++i) {
        array.push_back(std::rand());
    }
    
    // Sorts `array` through argsorting.
    std::vector<int> indices(array.size());
    argsort::argsort(indices.begin(), indices.end(), array.begin());
    std::vector<int> sorted_array(array.size());
    argsort::apply(
        indices.begin(), indices.end(), array.begin(), sorted_array.begin()
    );
    
    // Sorts `array` using std::sort and check it gives the same result.
    std::sort(array.begin(), array.end());
    assert(std::equal(array.begin(), array.end(), sorted_array.begin()));
}
