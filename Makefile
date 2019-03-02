CPP = c++ -std=c++11 -I/usr/local/include/eigen3

tests = test_argsort test_ConvexPolygon test_LinearFit test_LinearInterp test_NDIndex test_tuple_iterator

.PHONY: test
test: $(tests)

$(tests): test_%: test_%.cpp %.hpp
	$(CPP) -o $@ $<
	$@
