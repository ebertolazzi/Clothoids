SRCS = Sources/*.cc
CXX  = c++

all:	$(SRCS)
	$(CXX) -I./Headers -o test1 test1.cc $(SRCS)
	$(CXX) -I./Headers -o test2 test2.cc $(SRCS)
	$(CXX) -I./Headers -o test3 test3.cc $(SRCS)
	$(CXX) -I./Headers -o test4 test4.cc $(SRCS)
	$(CXX) -I./Headers -o test5 test5.cc $(SRCS)

run:
	./test1
	./test2
	./test3
	./test4
	./test5

doc:
	doxygen