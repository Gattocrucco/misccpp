#ifndef UNCERTAINTIES_IMPL_HPP_8AC76B25
#define UNCERTAINTIES_IMPL_HPP_8AC76B25

#include <string>

#include "ureal.hpp"

namespace uncertainties {
    Id next_id {};

    void _insert_dot(std::string *s, int n, int e) {
        e += s->size() - n;
        n = s->size();
        if (e >= n - 1) {
            // no dot at end of mantissa
        } else if (e >= 0) {
            s->insert(1 + e, 1, '.');
        } else if (e <= -1) {
            s->insert(0, -e, '0');
            s->insert(1, 1, '.');
        }
    }
}

#endif /* end of include guard: UNCERTAINTIES_IMPL_HPP_8AC76B25 */
