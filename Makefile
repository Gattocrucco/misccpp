CPP = c++
CFLAGS = -std=c++11 -g
LDFLAGS = -lc++

tests = test_argsort test_ConvexPolygon test_LinearFit test_LinearInterp test_NDIndex test_tuple_iterator

run_tests = $(patsubst test_%,run_%,$(tests))

.PHONY: test $(run_tests)
test: $(run_tests)

$(run_tests): run_%: test_%
	$<

$(tests): test_%: test_%.cpp %.hpp
	$(CPP) $(CFLAGS) $(LDFLAGS) -o $@ $<
test_LinearFit: CFLAGS += -I/usr/local/include/eigen3
