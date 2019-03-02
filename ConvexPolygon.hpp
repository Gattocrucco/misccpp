#ifndef CONVEXPOLYGON_HPP_96EBE97E
#define CONVEXPOLYGON_HPP_96EBE97E

#include <vector>
#include <stdexcept>
#include <cmath>

#include "tuple_iterator.hpp"

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
    
    ConvexPolygon(const std::vector<std::pair<Number, Number>> &vertices):
        ConvexPolygon(
            tuple_iterator::tuple_iterator<0>(vertices.begin()),
            tuple_iterator::tuple_iterator<0>(vertices.end()),
            tuple_iterator::tuple_iterator<1>(vertices.begin())
        )
    {;}
    
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
    
    std::vector<Number> copy_x() const {
        return get_from_vertices([](const Vertex &v) { return v.x; });
    }
    
    std::vector<Number> copy_y() const {
        return get_from_vertices([](const Vertex &v) { return v.y; });
    }
};

#endif /* end of include guard: CONVEXPOLYGON_HPP_96EBE97E */
