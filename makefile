SRCS = \
srcs/GenericContainer.cc \
srcs/GenericContainerSupport.cc \
srcs/GenericContainerCinterface.cc \
srcs_lua_interface/GenericContainerLuaInterface.cc

OBJS = $(SRCS:.cc=.o)
DEPS = \
srcs/GenericContainer.hh \
srcs/GenericContainerCinterface.h \
srcs_lua_interface/GenericContainerLuaInterface.hh

#CC     = llvm-gcc
#CXX    = llvm-g++
#CC     = clang
#CXX    = clang++
CC     = gcc
CXX    = g++

CFLAGS = -I./srcs -I./srcs_lua_interface -Wall -O3
LIBS   = -L./libs -lGenericContainer

#AR     = ar rcs
AR     = libtool -static -o 

all: libGenericContainer.a
	$(CXX) $(CFLAGS) -o bin/example1 examples/example1.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example2 examples/example2.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example3 examples/example3.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example4 examples/example4.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example5 examples/example5.cc $(LIBS)
	$(CC)  $(CFLAGS) -o bin/example6 examples/example6.c  $(LIBS) -lstdc++
	$(CXX) $(CFLAGS) -o bin/example7 examples/example7.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example8 examples/example8.cc $(LIBS) -llua
	$(CXX) $(CFLAGS) -o bin/example9 examples/example9.cc $(LIBS) -llua

srcs/%.o: srcs/%.cc $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@ 

srcs_lua_interface/%.o: srcs_lua_interface/%.cc $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@ 

srcs/%.o: srcs/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libGenericContainer.a: $(OBJS)
	$(AR) libs/libGenericContainer.* $(OBJS) 

run:
	cd bin ; ./example1
	cd bin ; ./example2
	cd bin ; ./example3
	cd bin ; ./example4
	cd bin ; ./example5
	cd bin ; ./example6
	cd bin ; ./example7
	cd bin ; ./example8
	cd bin ; ./example9

doc:
	doxygen
	
clean:
	rm -f libs/libGenericContainer.* srcs*/*.o
	rm -f bin/example1
	rm -f bin/example2
	rm -f bin/example3
	rm -f bin/example4
	rm -f bin/example5
	rm -f bin/example6
	rm -f bin/example7
	rm -f bin/example8
	rm -f bin/example9
	
