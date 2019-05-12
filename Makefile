CPP = c++
CFLAGS = -std=c++11 -g
LDFLAGS = -lc++

tests = test_argsort test_ConvexPolygon test_LinearFit test_LinearInterp test_NDIndex test_tuple_iterator
uncertainties_tests = test_uncertainties test_unc_urealhpp test_unc_iohpp test_unc_stathpp  test_unc_mathhpp test_unc_format test_unc_urealshpp test_unc_ureals

uncertainties_objects = $(patsubst test_%,test_%.cpp.o,$(uncertainties_tests))
run_tests = $(patsubst test_%,run_%,$(tests) $(uncertainties_tests))

.PHONY: test $(run_tests)
test: $(run_tests)

$(run_tests): run_%: test_%
	$<

$(tests): test_%: test_%.cpp %.hpp
	$(CPP) $(CFLAGS) $(LDFLAGS) -o $@ $<
test_LinearFit: CFLAGS += -I/usr/local/include/eigen3

$(uncertainties_tests): test_%: test_%.cpp.o uncertainties.cpp.o
	$(CPP) $(CFLAGS) $(LDFLAGS) -o $@ $^
uncertainties.cpp.o: uncertainties.cpp ./uncertainties/*.hpp
	$(CPP) $(CFLAGS) -I. -c -o $@ $<
$(uncertainties_objects): test_%.cpp.o: test_%.cpp ./uncertainties/*.hpp
	$(CPP) $(CFLAGS) -I. -c -o $@ $<
test_unc_urealshpp.cpp.o test_unc_ureals.cpp.o: CFLAGS += -I/usr/local/include/eigen3
