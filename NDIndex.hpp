#ifndef NDINDEX_HPP_C982B01F
#define NDINDEX_HPP_C982B01F

#include <cstddef>
#include <map>

// Multidimensional index iterator. See test_ndindex.cpp for an example.
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

    explicit NDIndex(std::map<key_t, size_t> sizes) {
        for (const auto key_size : sizes) {
            const key_t key = key_size.first;
            const size_t size = key_size.second;
            indices[key] = index_t{0, size};
        }
    }
    
    size_t operator[](key_t key) const {
        return indices.at(key).idx;
    }
    
    // Find next combination.
    // It is like adding 1 to a little endian number where each digit
    // uses a different base.
    NDIndex<key_t> &operator++() {
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
    
    operator bool() const {
        return not is_finished;
    }
};

#endif /* end of include guard: NDINDEX_HPP_C982B01F */
