# get the type of OS currently running
OS=$(shell uname -s)
PWD=$(shell pwd)

INC         = -I./src -I./include
LIBS        = -L./lib -lClothoids
DEFS        =
STATIC_EXT  = .a
DYNAMIC_EXT = .so
AR          = ar rcs
LDCONFIG    = sudo ldconfig

WARN=-Wall -Wno-sign-compare
#-Weverything -Wno-global-constructors -Wno-padded -Wno-documentation-unknown-command 

# default values

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  LIBS     = -static -L./lib -lClothoids
  CXXFLAGS = -std=c++11 $(WARN) -O3 -fPIC
  AR       = ar rcs
  LDCONFIG = sudo ldconfig
endif

# check if the OS string contains 'MINGW'
ifneq (,$(findstring MINGW, $(OS)))
  LIBS     = -static -L./lib -lClothoids
  CXXFLAGS = -std=c++11 $(WARN) -O3
  AR       = ar rcs
  LDCONFIG = sudo ldconfig
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  WARN        = -Wall -Weverything -Wno-sign-compare -Wno-global-constructors -Wno-padded -Wno-documentation-unknown-command 
  LIBS        = -L./lib -lClothoids
  CXXFLAGS    = $(WARN) -O3 -fPIC
  AR          = libtool -static -o
  LDCONFIG    =
  DYNAMIC_EXT = .dylib
endif

LIB_CLOTHOID = libClothoids

SRCS = \
src/AABBtree.cc \
src/Biarc.cc \
src/Circle.cc \
src/Clothoid.cc \
src/ClothoidDistance.cc \
src/ClothoidG2.cc \
src/ClothoidList.cc \
src/CubicRootsFlocke.cc \
src/Fresnel.cc \
src/G2lib.cc \
src/Line.cc \
src/Triangle2D.cc \
src/PolyLine.cc

OBJS  = $(SRCS:.cc=.o)
DEPS  = src/Clothoid.hh src/CubicRootsFlocke.hh
MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Clothoids

all: bin

travis: bin

bin: lib
	@$(MKDIR) bin
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testBiarc        tests-cpp/testBiarc.cc      $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testDistance     tests-cpp/testDistance.cc   $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2           tests-cpp/testG2.cc         $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2plot       tests-cpp/testG2plot.cc     $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2stat       tests-cpp/testG2stat.cc     $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2stat2arc   tests-cpp/testG2stat2arc.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2statCLC    tests-cpp/testG2statCLC.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testIntersect    tests-cpp/testIntersect.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testPolyline     tests-cpp/testPolyline.cc   $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testTriangle2D   tests-cpp/testTriangle2D.cc $(LIBS)

lib: lib/$(LIB_CLOTHOID)$(STATIC_EXT) lib/$(LIB_CLOTHOID)$(DYNAMIC_EXT)

include_local:
	@rm -rf lib/include
	$(MKDIR) lib
	$(MKDIR) lib/include
	@cp -f src/*.hh lib/include

src/%.o: src/%.cc $(DEPS)
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@ 

src/%.o: src/%.c $(DEPS)
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

lib/libClothoids.a: $(OBJS) include_local
	@$(MKDIR) lib
	$(AR) lib/libClothoids.a $(OBJS) 

lib/libClothoids.dylib: $(OBJS) include_local
	@$(MKDIR) lib
	$(CXX) -shared -o lib/libClothoids.dylib $(OBJS) 

lib/libClothoids.so: $(OBJS) include_local
	@$(MKDIR) lib
	$(CXX) -shared -o lib/libClothoids.so $(OBJS) 

install: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include
	cp src/*.hh                $(PREFIX)/include
	cp lib/$(LIB_CLOTHOID).*   $(PREFIX)/lib
	@$(LDCONFIG) $(PREFIX)/lib

install_as_framework: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp src/*.hh                $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_CLOTHOID)     $(PREFIX)/lib

run:
	./bin/testBiarc
	./bin/testDistance
	./bin/testG2
	./bin/testG2plot
	./bin/testG2stat
	./bin/testG2stat2arc
	./bin/testG2statCLC
	./bin/testIntersect
	./bin/testPolyline
	./bin/testTriangle2D

doc:
	doxygen
	
clean:
	rm -f lib/libClothoids.* lib/libClothoids.* src/*.o
	rm -rf bin
	
