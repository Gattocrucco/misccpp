# misccpp

Miscellaneous code in C++.

## Usage

They are just header files so place them alongside your code. Use `make` to run the tests.

## Documentation

None. See the tests for examples.

## Headers

* `argsort.hpp`: implements something like Numpy's `argsort`.

* `ConvexPolygon.hpp`: check if a polygon is convex and if a point is inside the polygon. Uses `tuple_iterator.hpp`.

* `LinearFit.hpp`: provides an interface to [least squares fitting](http://eigen.tuxfamily.org/dox/group__LeastSquares.html). Requires [Eigen](https://eigen.tuxfamily.org).

* `LinearInterp.hpp`: linear interpolator. Uses `argsort.hpp`.

* `NDIndex.hpp`: something like Numpy's `ndindex`, iterates over all combinations of indices.

* `tuple_iterator.hpp`: given an iterator of e.g. `std::tuple`, iterates over a fixed component of the tuples.
