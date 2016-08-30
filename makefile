# get the type of OS currently running
OS     = $(shell uname)
LIB_GC = libGenericContainer.a
CC     = gcc
CXX    = g++

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  #LIB_GC = libGenericContainer.so
  AR = ar rcs
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  CC  = clang   -Weverything -Wno-reserved-id-macro -Wno-padded
  CXX = clang++ -Weverything -Wno-reserved-id-macro -Wno-padded
  AR  = libtool -static -o
  #LIB_GC = libGenericContainer.dylib
endif

# to compile with lua make LUA_SUPPORT="YES"
LUALIB  =
SRCSLUA = 
ifneq (,$(findstring YES, $(LUA_SUPPORT)))
  SRCSLUA = \
  src_lua_interface/GenericContainerLuaInterface.cc \
  src_lua_interface/GenericContainerLuaPmain.cc
  LUALIB  = -llua
endif

SRCS = \
src/GenericContainer.cc \
src/GenericContainerSupport.cc \
src/GenericContainerTables.cc \
src/GenericContainerCinterface.cc \
$(SRCSLUA)

OBJS = $(SRCS:.cc=.o)
DEPS = \
src/GenericContainer.hh \
src/GenericContainerCinterface.h \
src_lua_interface/GenericContainerLuaInterface.hh

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = GenericContainer

CFLAGS   = -Wall -O3
CXXFLAGS = -Wall -O3
INC      = -I/usr/local/include -I./src -I./src_lua_interface
LIB_DIR  = -L/usr/local/lib -L./lib
LIBS     = $(LIB_DIR) -lGenericContainer -lpcre
DEFINE   =

MKDIR  = mkdir -p

all: lib
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example1  examples/example1.cc  $(LIBS)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example2  examples/example2.cc  $(LIBS)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example3  examples/example3.cc  $(LIBS)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example4  examples/example4.cc  $(LIBS)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example5  examples/example5.cc  $(LIBS)
	$(CC)  $(CFLAGS)   $(INC) -o bin/example6  examples/example6.c   $(LIBS) -lstdc++
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example7  examples/example7.cc  $(LIBS)
ifneq (,$(findstring lua, $(HASLUA)))
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example8  examples/example8.cc  $(LIBS) $(LUALIB)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example9  examples/example9.cc  $(LIBS) $(LUALIB)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example10 examples/example10.cc $(LIBS) $(LUALIB)
endif
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example11 examples/example11.cc $(LIBS)

lib: lib/$(LIB_GC)

src/%.o: src/%.cc $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@ 

src/%.o: src/%.c $(DEPS)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

src_lua_interface/%.o: src_lua_interface/%.cc $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) $(DEFINE) -c $< -o $@ 

lib/libGenericContainer.a: $(OBJS)
	$(MKDIR) lib
	$(AR) lib/libGenericContainer.a $(OBJS) 

lib/libGenericContainer.dylib: $(OBJS)
	$(MKDIR) lib
	$(CXX) -dynamiclib $(OBJS) -o lib/libGenericContainer.dylib $(LIB_DIR) -llua -lpcre -install_name libGenericContainer.dylib -Wl,-rpath,.

lib/libGenericContainer.so: $(OBJS)
	$(MKDIR) lib
	$(CXX) -shared $(OBJS) -o lib/libGenericContainer.so $(LIB_DIR) -llua -lpcre

install: lib/$(LIB_GC)
	cp src/GenericContainer.hh $(PREFIX)/include
	cp src_lua_interface/GenericContainerLuaInterface.hh $(PREFIX)/include
	cp lib/$(LIB_GC)           $(PREFIX)/lib

install_as_framework: lib/$(LIB_GC)
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp src/GenericContainer.hh $(PREFIX)/include/$(FRAMEWORK)
	cp src_lua_interface/GenericContainerLuaInterface.hh $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_GC)           $(PREFIX)/lib

run:
	cd bin ; ./example1
	cd bin ; ./example2
	cd bin ; ./example3
	cd bin ; ./example4
	cd bin ; ./example5
	cd bin ; ./example6
	cd bin ; ./example7
ifneq (,$(findstring lua, $(HASLUA)))
	cd bin ; ./example8
	cd bin ; ./example9
	cd bin ; ./example10
endif
	cd bin ; ./example11

doc:
	doxygen
	
clean:
	rm -f lib/libGenericContainer.* src*/*.o
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
