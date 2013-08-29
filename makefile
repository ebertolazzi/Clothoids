SRCS = \
Sources/GenericContainer.cc \
Sources/GenericContainerCinterface.cc \
Sources/GenericContainerSupport.cc

OBJS = $(SRCS:.cc=.o)
DEPS = Headers/GenericContainer.hh Headers/GenericContainerCinterface.h

#CC     = llvm-gcc
#CXX    = llvm-g++
#CC     = clang
#CXX    = clang++
CC     = gcc
CXX    = g++

CFLAGS =  -I./Headers -Wall -O3
LIBS   = -L. -lGenericContainer

#AR     = ar rcs
AR     = libtool -static -o 

all: libGenericContainer.a
	$(CXX) $(CFLAGS) -o test1 test1.cc $(LIBS)
	$(CXX) $(CFLAGS) -o test2 test2.cc $(LIBS)
	$(CXX) $(CFLAGS) -o test3 test3.cc $(LIBS)
	$(CXX) $(CFLAGS) -o test4 test4.cc $(LIBS)
	$(CXX) $(CFLAGS) -o test5 test5.cc $(LIBS)
	$(CC)  $(CFLAGS) -o test6 test6.c  $(LIBS) -lstdc++
	$(CXX) $(CFLAGS) -o test7 test7.cc $(LIBS)

Sources/%.o: Sources/%.cc $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@ 

Sources/%.o: Sources/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libGenericContainer.a: $(OBJS)
	$(AR) libGenericContainer.a $(OBJS) 

run:
	./test1
	./test2
	./test3
	./test4
	./test5
	./test6
	./test7

doc:
	doxygen
	
clean:
	rm -f libGenericContainer.a Sources/*.o