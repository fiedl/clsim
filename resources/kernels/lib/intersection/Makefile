all: test
run: test

test: clean test_cpu test_gpu

test_cpu: intersection_tests
	./intersection_tests

test_gpu: intersection_opencl_tests
	./intersection_opencl_tests

intersection_tests: intersection_tests.c intersection_tests.h intersection.c
	gcc \
		-Wno-implicit-function-declaration \
		-o intersection_tests \
		intersection_tests.c

intersection_opencl_tests: intersection_opencl_tests.c intersection_opencl_tests_kernel.cl intersection.c
	gcc \
	  -o intersection_opencl_tests \
		intersection_opencl_tests.c \
		-framework opencl

clean:
	touch intersection_tests; rm intersection_tests
	touch intersection_opencl_tests; rm intersection_opencl_tests
