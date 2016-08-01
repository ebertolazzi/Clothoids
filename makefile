# get the type of OS currently running
OS=$(shell uname)

LIB_GC = libGenericContainer.a
CC     = gcc
CXX    = g++

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  #LIB_GC = libGenericContainer.so
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  CC     = clang
  CXX    = clang++
  #LIB_GC = libGenericContainer.dylib
endif

SRCS = \
srcs/GenericContainer.cc \
srcs/GenericContainerSupport.cc \
srcs/GenericContainerTables.cc \
srcs/GenericContainerCinterface.cc \
srcs_lua_interface/GenericContainerLuaInterface.cc

OBJS = $(SRCS:.cc=.o)
DEPS = \
srcs/GenericContainer.hh \
srcs/GenericContainerCinterface.h \
srcs_lua_interface/GenericContainerLuaInterface.hh

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX=/usr/local

CFLAGS  = -I/usr/local/include -I./srcs -I./srcs_lua_interface -Wall -O3
LIB_DIR = -L/usr/local/lib -L./libs
LIBS    = $(LIB_DIR) -lGenericContainer -lpcre

#AR     = ar rcs
AR     = libtool -static -o
MKDIR  = mkdir -p

all: libs/$(LIB_GC)
	$(CXX) $(CFLAGS) -o bin/example1  examples/example1.cc  $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example2  examples/example2.cc  $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example3  examples/example3.cc  $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example4  examples/example4.cc  $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example5  examples/example5.cc  $(LIBS)
	$(CC)  $(CFLAGS) -o bin/example6  examples/example6.c   $(LIBS) -lstdc++
	$(CXX) $(CFLAGS) -o bin/example7  examples/example7.cc  $(LIBS)
	$(CXX) $(CFLAGS) -o bin/example8  examples/example8.cc  $(LIBS) -llua
	$(CXX) $(CFLAGS) -o bin/example9  examples/example9.cc  $(LIBS) -llua
	$(CXX) $(CFLAGS) -o bin/example10 examples/example10.cc $(LIBS) -llua
	$(CXX) $(CFLAGS) -o bin/example11 examples/example11.cc $(LIBS)

srcs/%.o: srcs/%.cc $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@ 

srcs_lua_interface/%.o: srcs_lua_interface/%.cc $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@ 

srcs/%.o: srcs/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libs/libGenericContainer.a: $(OBJS)
	$(MKDIR) libs
	$(AR) libs/libGenericContainer.a $(OBJS) 

libs/libGenericContainer.dylib: $(OBJS)
	$(MKDIR) libs
	$(CXX) -dynamiclib $(OBJS) -o libs/libGenericContainer.dylib $(LIB_DIR) -llua -lpcre -install_name libGenericContainer.dylib -Wl,-rpath,.

libs/libGenericContainer.so: $(OBJS)
	$(MKDIR) libs
	$(CXX) -shared $(OBJS) -o libs/libGenericContainer.so $(LIB_DIR) -llua -lpcre

install: libs/$(LIB_GC)
	cp srcs/GenericContainer.hh $(PREFIX)/include
	cp libs/$(LIB_GC)           $(PREFIX)/lib

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
	cd bin ; ./example10
	cd bin ; ./example11

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
	rm -f bin/example10
	rm -f bin/example11
