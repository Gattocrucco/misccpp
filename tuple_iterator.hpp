#ifndef TUPLE_ITERATOR_HPP_858B9520
#define TUPLE_ITERATOR_HPP_858B9520

namespace tuple_iterator {
    template<int n, typename Iterator>
    class TupleIterator {
    private:
        Iterator iterator;
        using Type = TupleIterator<n, Iterator>;
    public:
        TupleIterator(Iterator it): iterator {it} {
            ;
        }
        decltype(std::get<n>(*iterator)) operator*() const {
            return std::get<n>(*iterator);
            // Note: `decltype` preserves references, and `std::get` has the
            // correct reference type based on `*iterator`.
        }
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

    template<int n, typename Iterator>
    TupleIterator<n, Iterator> tuple_iterator(Iterator it) {
        return it;
    }
}

#endif /* end of include guard: TUPLE_ITERATOR_HPP_858B9520 */
