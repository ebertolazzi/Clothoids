# get the type of OS currently running
OS=$(shell uname -s)
PWD=$(shell pwd)

INC         = -I./src -I./include -Isubmodules/quarticRootsFlocke/src
DEFS        =
STATIC_EXT  = .a
DYNAMIC_EXT = .so
AR          = ar rcs

WARN        = -Wall -Wno-sign-compare
#-Weverything -Wno-global-constructors -Wno-padded -Wno-documentation-unknown-command

# default values
LIB_CLOTHOID = Clothoids
LIBS         = -L./lib/lib -l$(LIB_CLOTHOID)_static
CXXFLAGS     = -O2
AR           = ar rcs
LDCONFIG     =


# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  LIB_CLOTHOID = Clothoids_linux
  LIBS         = -L./lib/lib -l$(LIB_CLOTHOID)_static
  CXXFLAGS     = -std=c++11 $(WARN) -O2 -fPIC
  AR           = ar rcs
  LDCONFIG     = sudo ldconfig
endif

# check if the OS string contains 'MINGW'
ifneq (,$(findstring MINGW, $(OS)))
  LIB_CLOTHOID = Clothoids_mingw_x64
  LIBS         = -L./lib/lib -l$(LIB_CLOTHOID)_static
  CXXFLAGS     = -std=c++11 $(WARN) -O2
  AR           = ar rcs
  LDCONFIG     = sudo ldconfig
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  LIB_CLOTHOID = Clothoids_osx
  LIBS         = -L./lib/lib -l$(LIB_CLOTHOID)_static
  WARN         = -Wall -Weverything -Wno-sign-compare -Wno-global-constructors -Wno-padded -Wno-documentation-unknown-command
	CC           = clang
	CXX          = clang++ -std=c++11
  CXXFLAGS     = $(WARN) -O2 -fPIC
  AR           = libtool -static -o
  LDCONFIG     =
  DYNAMIC_EXT  = .dylib
endif

.SUFFIXES: .o

SRCS = \
src/AABBtree.cc \
src/Biarc.cc \
src/BiarcList.cc \
src/Circle.cc \
src/Clothoid.cc \
src/ClothoidAsyPlot.cc \
src/ClothoidDistance.cc \
src/ClothoidG2.cc \
src/ClothoidList.cc \
src/Fresnel.cc \
src/G2lib.cc \
src/G2lib_intersect.cc \
src/Line.cc \
src/PolyLine.cc \
src/Triangle2D.cc \
submodules/quarticRootsFlocke/src/PolynomialRoots-1-Quadratic.cc \
submodules/quarticRootsFlocke/src/PolynomialRoots-2-Cubic.cc \
submodules/quarticRootsFlocke/src/PolynomialRoots-3-Quartic.cc \
submodules/quarticRootsFlocke/src/PolynomialRoots-Utils.cc

OBJS  = $(SRCS:.cc=.o)
DEPS  = src/Clothoid.hh src/CubicRootsFlocke.hh
MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Clothoids

all: bin

travis: bin run

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

lib: lib/lib/lib$(LIB_CLOTHOID)_static$(STATIC_EXT) lib/lib/lib$(LIB_CLOTHOID)$(DYNAMIC_EXT)

include_local:
	@rm -rf lib/include
	@$(MKDIR) -p lib/include
	@cp -f src/*.hh                               lib/include
	@cp -f submodules/quarticRootsFlocke/src/*.hh lib/include


.cc.o : $(DEPS)
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

.c.o : $(DEPS)
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

lib/lib/lib$(LIB_CLOTHOID)_static.a: $(OBJS) include_local
	@$(MKDIR) -p lib/lib
	$(AR) lib/lib/lib$(LIB_CLOTHOID)_static.a $(OBJS)

lib/lib/lib$(LIB_CLOTHOID).dylib: $(OBJS) include_local
	@$(MKDIR) -p lib/lib
	$(CXX) -shared -o lib/lib/lib$(LIB_CLOTHOID).dylib $(OBJS)

lib/lib/lib$(LIB_CLOTHOID).so: $(OBJS) include_local
	@$(MKDIR) -p lib/lib
	$(CXX) -shared -o lib/lib/lib$(LIB_CLOTHOID).so $(OBJS)

install: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include
	@cp src/*.hh                                   $(PREFIX)/include
	@cp src/submodules/quarticRootsFlocke/src/*.hh $(PREFIX)/include
	@cp lib/lib/lib$(LIB_CLOTHOID).*               $(PREFIX)/lib
	@$(LDCONFIG) $(PREFIX)/lib

install_as_framework: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	@cp src/*.hh                                   $(PREFIX)/include/$(FRAMEWORK)
	@cp src/submodules/quarticRootsFlocke/src/*.hh $(PREFIX)/include/$(FRAMEWORK)
	@cp lib/lib/lib$(LIB_CLOTHOID)                 $(PREFIX)/lib

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

docs:
	@doxygen
	@open docs/index.html

clean:
	rm -f lib/libClothoids.* src/*.o
	rm -rf bin
