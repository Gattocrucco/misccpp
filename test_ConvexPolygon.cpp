#include <vector>
#include <stdexcept>
#include <iostream>
#include <cassert>

#include "ConvexPolygon.hpp"

using std::vector;
using std::exception;
using std::cerr;
using std::pair;

template<typename Operation>
bool throws(Operation operation) {
    bool throwed;
    try {
        operation();
        throwed = false;
    } catch(exception &e) {
        throwed = true;
    } catch(...) {
        cerr << "throws(): cannot catch exception";
        throw;
    }
    return throwed;
}

int main() {
    // arrow
    assert(throws([]() {
        ConvexPolygon<double> polygon({
            {0, 0},
            {1, 1},
            {-1, 0},
            {1, -1}
        });
    }));
    
    // glasshour
    assert(throws([]() {
        ConvexPolygon<double> polygon({
            {0, 0},
            {1, 1},
            {1, 0},
            {0, 1},
        });
    }));
    
    // boomerang
    assert(throws([]() {
        ConvexPolygon<double> polygon({
            {0, 0},
            {3, 0},
            {3, 2},
            {2, 1}
        });
    }));
    
    // twisted pentagon
    assert(throws([]() {
        ConvexPolygon<double> polygon({
            {0, 0},
            {4, 2},
            {2, 2},
            {6, 0},
            {3, 4},
        });
    }));
    
    // clockwise square
    ConvexPolygon<double> p1({
        {0, 0},
        {0, 1},
        {1, 1},
        {1, 0}
    });
    // counter-clockwise square
    ConvexPolygon<double> p2({
        {0, 0},
        {1, 0},
        {1, 1},
        {0, 1}
    });
    vector<pair<double, double>> points_inside {
        {0.5, 0.5},
        {0.1, 0.1},
        {0.1, 0.9},
        {0, 0},
        {0, 1},
        {1, 1},
        {1, 0},
        {0.5, 0},
        {1, 0.5},
        {0.5, 1},
        {0, 0.5}
    };
    vector<pair<double, double>> points_outside {
        {1.5, 0.5},
        {0.5, 1.5},
        {1.5, 1.5},
        {0, 1.5},
        {1, 1.5},
        {-0.5, 1},
        {0, -0.5},
    };
    assert(std::all_of(points_inside.begin(), points_inside.end(),
    [&](pair<double, double> p) {
        return p1.contains_point(p.first, p.second);
    }));
    assert(std::all_of(points_inside.begin(), points_inside.end(),
    [&](pair<double, double> p) {
        return p2.contains_point(p.first, p.second);
    }));
    assert(not std::any_of(points_outside.begin(), points_outside.end(),
    [&](pair<double, double> p) {
        return p1.contains_point(p.first, p.second);
    }));
    assert(not std::any_of(points_outside.begin(), points_outside.end(),
    [&](pair<double, double> p) {
        return p2.contains_point(p.first, p.second);
    }));
}
