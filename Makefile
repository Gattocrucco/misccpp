CPP = c++ -std=c++11 -I/usr/local/include/eigen3 -g

tests = test_argsort test_ConvexPolygon test_LinearFit test_LinearInterp test_NDIndex test_tuple_iterator test_uncertainties
run_tests = $(patsubst test_%,run_%,$(tests))

.PHONY: test $(run_tests)
test: $(run_tests)

$(run_tests): run_%: test_%
	$<

$(tests): test_%: test_%.cpp %.hpp
	$(CPP) -o $@ $<
