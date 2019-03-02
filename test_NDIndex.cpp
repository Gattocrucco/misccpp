#include <iostream>
#include <vector>
#include <cstddef>
#include <map>

#include "NDIndex.hpp"

using namespace std;

int main() {
    vector<float> v1 {1, 2};
    vector<float> v2 {3, 4, 5};
    vector<float> v3 {6, 7, 8, 9};
    
    map<int, size_t> sizes {
        {1, v1.size()},
        {2, v2.size()},
        {3, v3.size()}
    };
    
    for (NDIndex<int> i(sizes); i; ++i) {
        cout << v1[i[1]] << ", " << v2[i[2]] << ", " << v3[i[3]] << "\n";
    }
}
