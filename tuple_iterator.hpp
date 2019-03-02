#ifndef TUPLE_ITERATOR_HPP_858B9520
#define TUPLE_ITERATOR_HPP_858B9520

/*! \file
\brief Defines class `TupleIterator` to iterate over a sequence of tuples
getting only a fixed element of each tuple.

Since `TupleIterator` uses `std::get`, you have to include the appropriate
headers before `tuple_iterator.hpp` (tipically `<tuple>`).

Example:
~~~{.cpp}
#include <iostream>
#include <vector>
#include <pair> // included before `tuple_iterator.hpp`
#include "tuple_iterator.hpp"
int main() {
    std::vector<std::pair<int, float>> v {
        {10, 3.14},
        {20, 6.28}
    };
    // it works for anything that `std::get` understands, so also `std::pair`s
    auto begin = tuple_iterator::tuple_iterator<0>(v.begin());
    auto end = tuple_iterator::tuple_iterator<0>(v.end());
    // we used the convenience function `tuple_iterator` to make
    // the `TupleIterator`s
    for (it = begin; it != end; ++it) {
        std::cout << *it << '\n'; // will print `10` then `20`
    }
}
~~~
*/

/*!
\brief Namespace for all the classes and functions defined in
`tuple_iterator.hpp`.
*/
namespace tuple_iterator {
    /*!
    \brief Class template that implements the tuple iterator.
    
    The template parameter `n` is the index of the tuple to get the element.
    The parameter `Iterator` is the type of the internal iterator to tuples,
    to avoid specifying it explicitly use the function
    `tuple_iterator::tuple_iterator`.
    */
    template<int n, typename Iterator>
    class TupleIterator {
    private:
        Iterator iterator;
        using Type = TupleIterator<n, Iterator>;
    public:
        /*!
        \brief Construct the iterator given an iterator of tuples.
        
        The iterated type does not really need to be a tuple, anything that
        can be `std::get`ted is fine.
        */
        TupleIterator(Iterator it): iterator {it} {
            ;
        }
        
        /*!
        \brief Dereference the iterator passed at construction and applies
        `std::get<n>` to it, where `n` is the class template parameter.
        */
        decltype(std::get<n>(*iterator)) operator*() const {
            return std::get<n>(*iterator);
            // Note: `decltype` preserves references, and `std::get` has the
            // correct reference type based on `*iterator`.
        }
        
        /*!
        \brief Increments the iterator.
        */
        Type &operator++() {
            ++iterator;
            return *this;
        }
        friend bool operator==(const Type &i1, const Type &i2) {
            return i1.iterator == i2.iterator;
        }
        friend bool operator!=(const Type &i1, const Type &i2) {
            return i1.iterator != i2.iterator;
        }
    };
    
    /*!
    \brief Convenience function to construct a `TupleIterator`.
    
    Using this function you only need to explicitly specify the tuple index `n`
    but not the type of the iterator of tuples.
    */
    template<int n, typename Iterator>
    TupleIterator<n, Iterator> tuple_iterator(Iterator it) {
        return it;
    }
}

#endif /* end of include guard: TUPLE_ITERATOR_HPP_858B9520 */
