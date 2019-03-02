#ifndef NDINDEX_HPP_C982B01F
#define NDINDEX_HPP_C982B01F

#include <cstddef>
#include <map>

/*! \file
\brief Defines class `NDIndex` to iterate over multidimensional arrays.
*/

/*!
\brief Class template to iterate over multidimensional arrays.

This is similar to Numpy's `ndindex`. It is useful to iterate over all the
elements of a multidimensional array or over all combinations of elements of
various arrays.

To be more flexible it uses a map to identify the various indices instead
of the standard `0...(d-1)`. This may be useful in the "all the combinations"
case. The template parameter `key_t`, saved as member type `key_type`, is the
type of the keys. Read the following example to understand this map thing:
~~~{.cpp}
#include <iostream>
#include "NDIndex.hpp"
int main() {
    int a[3][4];
    NDIndex<int> idx(std::map{{10, 3}, {20, 4}});
    // the map given at construction says that the sizes are 3 and 4; the keys
    // `10` and `20` are arbitrary
    for (; idx; ++idx) {
        // the index will evaluate to `false` when it is finished
        std::cout << a[idx[10]][idx[20]] << '\n'; // prints the contents of `a`
        std::cout << idx[10] << ', ' << idx[20] << '\n'; // will print:
        // 0, 0
        // 0, 1
        // ...
        // 0, 3
        // 1, 0
        // 1, 1
        // ... etc
        // 2, 3
    }
}
~~~
*/
template<typename key_t>
class NDIndex {
private:    
    struct index_t {
        size_t idx;
        size_t size;
    };
    
    // IMPORTANT: use `map`, not `unordered_map`.
    std::map<key_t, index_t> indices;
    bool is_finished {false};

public:
    using key_type = key_t;
    
    /*!
    \brief Constructs the index given a map telling the sizes of the indices.
    
    The choice of the keys of the map is totally up to the user; they will be
    needed to extract the individual indices from the index using `operator[]`.
    Each individual index will iterate up to to its size - 1.
    */
    explicit NDIndex(std::map<key_t, size_t> sizes) {
        for (const auto key_size : sizes) {
            const key_t key = key_size.first;
            const size_t size = key_size.second;
            indices[key] = index_t{0, size};
        }
    }
    
    /*!
    \brief Get the value of one of the individual indices.
    
    The key must be one of the keys in the map given at construction.
    */
    size_t operator[](key_t key) const {
        return indices.at(key).idx;
    }
    
    /*!
    \brief Increments the index.
    
    The individual indices are incremented lexicographically. The order of the
    indices is given by the order of the keys (recall that maps are ordered).
    */
    NDIndex<key_t> &operator++() {
        // It is like adding 1 to a little endian number where each digit
        // uses a different base.
        size_t zero_count = 0;
        for (auto &key_index : indices) {
            // IMPORTANT: use references because we want to modify objects
            // inside `indices`.
            size_t &idx = key_index.second.idx;
            size_t &size = key_index.second.size;
            ++idx;
            if (idx >= size) {
                idx = 0;
                ++zero_count;
            } else {
                break;
            }
        }
        if (zero_count == indices.size()) {
            is_finished = true;
        }
        return *this;
    }
    
    /*!
    \brief Automatic cast to bool, evaluates to `false` when all the possible
    combinations of indices have been iterated over.
    */
    operator bool() const {
        return not is_finished;
    }
};

#endif /* end of include guard: NDINDEX_HPP_C982B01F */
