GTEST_DIR = ../gtest
USER_DIR = .
CPPFLAGS += -isystem $(GTEST_DIR)/include
CXXFLAGS += -g -Wall -Wextra -pthread
TESTS = intersection_test
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
								$(GTEST_DIR)/include/gtest/internal/*.h

all : $(TESTS) test
clean :
	rm -f $(TESTS) gtest.a gtest_main.a *.o

GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
						$(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
						$(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
		$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
		$(AR) $(ARFLAGS) $@ $^

#intersection.o : $(USER_DIR)/intersection.c $(USER_DIR)/intersection.h $(GTEST_HEADERS)
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/intersection.c

intersection_test.o : $(USER_DIR)/intersection_test.c \
										 $(USER_DIR)/intersection.c $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/intersection_test.c

intersection_test : intersection_test.o gtest_main.a
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@

intersection_test_opencl.o : $(USER_DIR)/intersection_test_opencl.c $(USER_DIR)/intersection.c $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/intersection_test_opencl.c

intersection_test_opencl: intersection_test_opencl.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ -framework opencl

test : intersection_test intersection_test_opencl
	$(USER_DIR)/intersection_test
	$(USER_DIR)/intersection_test_opencl
