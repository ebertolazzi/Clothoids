# get the type of OS currently running
OS      = $(shell uname)
LIB_GC  = libGenericContainer.a
CC      = gcc
CXX     = g++
INC    += -I./lib3rd/include
LIBSGCC =

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  WARN = -Wall
  CC  = gcc
  CXX = g++
  # activate C++11 for g++ >= 4.9
  VERSION  = $(shell $(CC) -dumpversion)
ifneq (,$(findstring 4.9, $(VERSION)))
  CXX += -std=c++11
endif
ifneq (,$(findstring 5., $(VERSION)))
  CXX += -std=c++11
endif
ifneq (,$(findstring 6., $(VERSION)))
  CXX += -std=c++11
endif
  CXXFLAGS = -pthread -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -m64 -O3 -g0 -funroll-loops -fPIC
  CFLAGS   = -pthread -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -m64 -O3 -g0 -funroll-loops -fPIC
  CC      += $(WARN)
  CXX     += $(WARN)
  AR       = ar rcs
  LIBSGCC = -lstdc++ -lm
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  WARN    = -Weverything -Wno-reserved-id-macro -Wno-padded
  CC      = clang
  CXX     = clang++
  VERSION = $(shell $(CC) --version 2>&1 | grep -o "Apple LLVM version [0-9]\.[0-9]\.[0-9]" | grep -o " [0-9]\.")
ifneq (,$(findstring 8., $(VERSION)))
  CXX += -std=c++11 -stdlib=libc++ 
endif
ifneq (,$(findstring 7., $(VERSION)))
  CXX += -std=c++11 -stdlib=libc++ 
endif
  CC      += $(WARN)
  CXX     += $(WARN)
  CXXFLAGS = -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -m64 -O3 -g0 -funroll-loops -fPIC
  CXXFLAGS = -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -m64 -O3 -g0 -funroll-loops -fPIC
  AR       = libtool -static -o
  LIBSGCC  = -lstdc++ -lm
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
  INC    += -Isrc_lua_interface
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

INC     += -I/usr/local/include -I./src -I./src_lua_interface
LIB_DIR  = -L/usr/local/lib -L./lib
LIBS     = $(LIB_DIR) -lGenericContainer -lpcre -lm
DEFINE   =

MKDIR  = mkdir -p

all: lib
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example1  examples/example1.cc  $(LIBS)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example2  examples/example2.cc  $(LIBS)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example3  examples/example3.cc  $(LIBS)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example4  examples/example4.cc  $(LIBS)
	$(CXX) $(CXXFLAGS) $(INC) -o bin/example5  examples/example5.cc  $(LIBS)
	$(CC)  $(CFLAGS)   $(INC) -o bin/example6  examples/example6.c   $(LIBS) $(LIBSGCC)
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

include_local:
	@rm -rf lib/include
	$(MKDIR) lib
	$(MKDIR) lib/include
	@cp -f src/GenericContainer.hh lib/include
	@cp -f src/GenericContainerConfig.hh lib/include
	@cp -f src/GenericContainerCinterface.h lib/include
	@cp -f src_lua_interface/GenericContainerLuaInterface.hh lib/include

lib/libGenericContainer.a: $(OBJS) include_local
	$(MKDIR) lib
	$(AR) lib/libGenericContainer.a $(OBJS) 

lib/libGenericContainer.dylib: $(OBJS) include_local
	$(MKDIR) lib
	$(CXX) -dynamiclib $(OBJS) -o lib/libGenericContainer.dylib $(LIB_DIR) -llua -lpcre -install_name libGenericContainer.dylib -Wl,-rpath,.

lib/libGenericContainer.so: $(OBJS) include_local
	$(MKDIR) lib
	$(CXX) -shared $(OBJS) -o lib/libGenericContainer.so $(LIB_DIR) -llua -lpcre

install: lib/$(LIB_GC)
	cp lib/include/GenericContainer.hh $(PREFIX)/include
	cp lib/include/GenericContainerConfig.hh $(PREFIX)/include
	cp lib/include/GenericContainerCinterface.h $(PREFIX)/include
	cp lib/include/GenericContainerLuaInterface.hh $(PREFIX)/include
	cp lib/$(LIB_GC) $(PREFIX)/lib

install_as_framework: lib/$(LIB_GC)
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp lib/include/GenericContainerConfig.hh $(PREFIX)/include/$(FRAMEWORK)
	cp lib/include/GenericContainerCinterface.h $(PREFIX)/include/$(FRAMEWORK)
	cp lib/include/GenericContainer.hh $(PREFIX)/include/$(FRAMEWORK)
	cp lib/include/GenericContainerLuaInterface.hh $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_GC) $(PREFIX)/lib

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
