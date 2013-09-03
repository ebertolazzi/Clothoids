SRCS = \
srcs/GenericContainer.cc \
srcs/GenericContainerCinterface.cc \
srcs/GenericContainerSupport.cc

OBJS = $(SRCS:.cc=.o)
DEPS = srcs/GenericContainer.hh srcs/GenericContainerCinterface.h

#CC     = llvm-gcc
#CXX    = llvm-g++
#CC     = clang
#CXX    = clang++
CC     = gcc
CXX    = g++

CFLAGS =  -I./srcs -Wall -O3
LIBS   = -L./libs -lGenericContainer

#AR     = ar rcs
AR     = libtool -static -o 

all: libGenericContainer.a
	$(CXX) $(CFLAGS) -o examples/example1 examples/example1.cc $(LIBS)
	$(CXX) $(CFLAGS) -o examples/example2 examples/example2.cc $(LIBS)
	$(CXX) $(CFLAGS) -o examples/example3 examples/example3.cc $(LIBS)
	$(CXX) $(CFLAGS) -o examples/example4 examples/example4.cc $(LIBS)
	$(CXX) $(CFLAGS) -o examples/example5 examples/example5.cc $(LIBS)
	$(CC)  $(CFLAGS) -o examples/example6 examples/example6.c  $(LIBS) -lstdc++
	$(CXX) $(CFLAGS) -o examples/example7 examples/example7.cc $(LIBS)

Sources/%.o: Sources/%.cc $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@ 

Sources/%.o: Sources/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libGenericContainer.a: $(OBJS)
	$(AR) libs/libGenericContainer.a $(OBJS) 

run:
	cd examples ; ./example1
	cd examples ; ./example2
	cd examples ; ./example3
	cd examples ; ./example4
	cd examples ; ./example5
	cd examples ; ./example6
	cd examples ; ./example7

doc:
	doxygen
	
clean:
	rm -f libs/libGenericContainer.a srcs/*.o
	rm -f examples/example1
	rm -f examples/example2
	rm -f examples/example3
	rm -f examples/example4
	rm -f examples/example5
	rm -f examples/example6
	rm -f examples/example7
	