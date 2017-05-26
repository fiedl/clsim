all: test
run: test

test: hole_ice_tests
	./hole_ice_tests

hole_ice_tests: hole_ice_tests.c hole_ice_tests.h hole_ice.c
	gcc \
		-Wno-implicit-function-declaration \
		-o hole_ice_tests \
		hole_ice_tests.c

clean:
	touch hole_ice_tests; rm hole_ice_tests
