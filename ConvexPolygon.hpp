#ifndef CONVEXPOLYGON_HPP_96EBE97E
#define CONVEXPOLYGON_HPP_96EBE97E

#include <vector>
#include <stdexcept>
#include <cmath>

#include "tuple_iterator.hpp"

/*! \file
\brief Defines class template `ConvexPolygon`.
*/

/*!
\brief Implements a convex polygon that can check if it contains a point.

A polygon can be checked to be convex in linear time, and if it is convex it
can be checked if it contains a given point in linear time (linear in the
number of vertices).

`ConvexPolygon` checks that the polygon is convex at construction. If it is not
it throws an exception.

The numerical type used to save the coordinates of the vertices and to make
computations is a template parameter. Both checking if the polygon is convex and
if it contains a point are exact computations if the numerical type is integer
(also if the type is floating-point but one uses only integer coordinates low
enough to avoid rounding). In the computations the numbers are multiplied so
one must be careful about overflows, which are not checked. Consider that each
coordinate is squared at some point.

Example:
~~~{.cpp}
#include "ConvexPolygon.hpp"
int main() {
    ConvexPolygon<double> polygon({
        // x, y
        {0, 0},
        {1, 0},
        {0.5, 1}
    }); // a triangle
    polygon.contains_point(0, 0); // returns `true` (edges count as inside)
    polygon.contains_point(2, 2); // returns `false`
}
~~~
*/ 
template<typename Number>
class ConvexPolygon {
private:
    struct Vertex {
        const Number x, y; // coordinates of vertex
        const Number vx, vy; // vector (arrow) to next vertex
        Vertex(
            const Number &this_x, const Number &this_y,
            const Number &next_x, const Number &next_y
        ):
            x {this_x}, y {this_y},
            vx {next_x - this_x}, vy {next_y - this_y}
        {;}
    };
    std::vector<Vertex> vertices;
    
    template<typename XIterator, typename YIterator>
    void fill_vertices(XIterator xbegin, XIterator xend, YIterator ybegin) {
        while (xbegin != xend) {
            const Number &x {*xbegin};
            const Number &y {*ybegin};
            ++xbegin;
            ++ybegin;
            if (xbegin == xend) {
                this->vertices.emplace_back(
                    x, y, vertices.front().x, vertices.front().y
                );
            } else {
                this->vertices.emplace_back(x, y, *xbegin, *ybegin);
            }
        }
    }
    
    static Number vector_product(const Vertex &a, const Vertex &b) {
        return vector_product(a.vx, a.vy, b.vx, b.vy);
    }
    
    static Number vector_product(
        const Number &ax, const Number &ay,
        const Number &bx, const Number &by
    ) {
        return ax * by - ay * bx;
    }
    
    static int sign(const Number &a) {
        if (a > 0) {
            return 1;
        } else if (a < 0) {
            return -1;
        } else {
            return 0;
        }
    }
    
    static bool same_sign_weak(
        const Number &a, const Number &b, const Number &c
    ) {
        return (a >= 0 and b >= 0 and c >= 0) or (a <= 0 and b <= 0 and c <= 0);
    }
    
    template<typename T>
    static bool different_sign_weak_right(const T &a, const T &b) {
        return (a > 0 and b <= 0) or (a < 0 and b >= 0);
    }
    
    bool is_convex() const {
        Number first_nonzero_prod = 0;
        size_t inversion_x_count = 0, inversion_y_count = 0;
        for (size_t i = 0; i < vertices.size(); ++i) {
            const Vertex &a = i > 0 ? vertices[i - 1] : vertices.back();
            const Vertex &b = vertices[i];
            const Number prod = vector_product(a, b);
            if (first_nonzero_prod != 0) {
                if (different_sign_weak_right(first_nonzero_prod, prod)) {
                    return false;
                }
            } else if (prod != 0) {
                first_nonzero_prod = std::move(prod);
            }
            if (different_sign_weak_right(a.vx, b.vx)) ++inversion_x_count;
            if (different_sign_weak_right(a.vy, b.vy)) ++inversion_y_count;
            if (inversion_x_count > 2 or inversion_y_count > 2) {
                return false;
            }
        }
        return true;
    }
    
    template<typename Getter>
    std::vector<Number> get_from_vertices(Getter getter) const {
        std::vector<Number> v;
        v.reserve(vertices.size());
        for (const Vertex &vertex : vertices) {
            v.push_back(getter(vertex));
        }
        return v;
    }
        
public:
    using coord_type = Number;
    
    /*!
    \brief Constructs the polygon from the coordinates of the vertices, given
    as two separate sequences for x and y.
    
    If the polygon is not convex, it throws `std::runtime_error`.
    */
    template<typename XIterator, typename YIterator>
    explicit ConvexPolygon(XIterator xbegin, XIterator xend, YIterator ybegin) {
        this->fill_vertices(xbegin, xend, ybegin);
        if (this->vertices.size() < 3) {
            throw std::runtime_error("ConvexPolygon: less than 3 vertices");
        }
        if (not this->is_convex()) {
            throw std::runtime_error(
                "ConvexPolygon: vertices do not form a convex polygon"
            );
        }
    }
    
    /*!
    \brief Constructs the polygon from the coordinates of the vertices, given
    as a vector of pairs {x, y}.
    
    If the polygon is not convex, it throws `std::runtime_error`.
    */
    ConvexPolygon(const std::vector<std::pair<Number, Number>> &vertices):
        ConvexPolygon(
            tuple_iterator::tuple_iterator<0>(vertices.begin()),
            tuple_iterator::tuple_iterator<0>(vertices.end()),
            tuple_iterator::tuple_iterator<1>(vertices.begin())
        )
    {;}
    
    /*!
    \brief Checks if the polygon contains a point.
    
    The point is specified by its coordinates `x` and `y`. Points on the edges
    of the polygon are considered inside. If `Number` (aka `coord_type`) is an
    integer type, the computation is exact (it does not suffer from rounding
    errors).
    */
    bool contains_point(const Number &x, const Number &y) const {
        for (size_t i = 0; i < vertices.size(); ++i) {
            const Vertex &curr = vertices[i];
            const Vertex &prev = i > 0 ? vertices[i - 1] : vertices.back();
            const Number vx {x - curr.x};
            const Number vy {y - curr.y};
            const Number prod_prevcurr = vector_product(prev, curr);
            const Number prod_prev = vector_product(prev.vx, prev.vy, vx, vy);
            const Number prod_curr = vector_product(curr.vx, curr.vy, vx, vy);
            if (not same_sign_weak(prod_prevcurr, prod_prev, prod_curr)) {
                return false;
            }
        }
        return true;
    }
    
    /*!
    \brief Return a vector of the x coordinate of the vertices, with the
    order given at construction.
    */
    std::vector<Number> copy_x() const {
        return get_from_vertices([](const Vertex &v) { return v.x; });
    }
    
    /*!
    \brief Return a vector of the y coordinate of the vertices, with the
    order given at construction.
    */
    std::vector<Number> copy_y() const {
        return get_from_vertices([](const Vertex &v) { return v.y; });
    }
};

#endif /* end of include guard: CONVEXPOLYGON_HPP_96EBE97E */
