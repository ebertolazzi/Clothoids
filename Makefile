# get the type of OS currently running
OS=$(shell uname -s)
PWD=$(shell pwd)

INC         = -I./src -I./include -Isubmodules/Utils/src -Isubmodules/Utils/src/Utils -Isubmodules/quarticRootsFlocke/src -Isubmodules/GenericContainer/src
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
  LIBS         = -L./lib/lib -l$(LIB_CLOTHOID)_static -ldl
  CXXFLAGS     = -std=c++11 $(WARN) -O2 -fPIC
  AR           = ar rcs
  LDCONFIG     = sudo ldconfig
endif

# check if the OS string contains 'MINGW'
ifneq (,$(findstring MINGW, $(OS)))
  LIB_CLOTHOID = Clothoids_mingw_x64
  LIBS         = -L./lib/lib -l$(LIB_CLOTHOID)_static
  CXXFLAGS     = -std=c++11 $(WARN) -O2 -Wsuggest-override
  AR           = ar rcs
  LDCONFIG     = sudo ldconfig
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  LIB_CLOTHOID = Clothoids_osx
  LIBS         = -L./lib/lib -l$(LIB_CLOTHOID)_static
  WARN         = -Wall -Wno-sign-compare -Wno-global-constructors -Wno-padded -Wno-documentation-unknown-command -Wno-poison-system-directories
	CC           = clang
	CXX          = clang++ -std=c++11
  CXXFLAGS     = $(WARN) -O2 -fPIC
  AR           = libtool -static -o
  LDCONFIG     =
  DYNAMIC_EXT  = .dylib
endif

.SUFFIXES: .o

SRCS  = $(shell echo src/*.cc) \
        $(shell echo submodules/Utils/src/*.cc) \
        $(shell echo submodules/quarticRootsFlocke/src/*.cc) \
        $(shell echo submodules/GenericContainer/src/*.cc)

OBJS  = $(SRCS:.cc=.o)
MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Clothoids

all: bin

travis: bin run

bin: lib
	@$(MKDIR) -p bin
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

lib: include_local lib/lib/lib$(LIB_CLOTHOID)_static$(STATIC_EXT) lib/lib/lib$(LIB_CLOTHOID)$(DYNAMIC_EXT)

include_local: .FORCE
	@rm -rf lib/include
	@$(MKDIR) -p lib/include/Utils
	@$(MKDIR) -p lib/include/Clothoids
	@cp -f src/*.hh                               lib/include
	@cp -f src/Clothoids/*.hxx                    lib/include/Clothoids

include_local_all: include_local
	@cp -f submodules/quarticRootsFlocke/src/*.hh lib/include/$(FRAMEWORK)
	@cp -f submodules/Utils/src/*.hh              lib/include/$(FRAMEWORK)
	@cp -rf submodules/Utils/src/Utils            lib/include/$(FRAMEWORK)
	@cp submodules/Utils/lib/lib/*.*              lib/lib
	@cp submodules/quarticRootsFlocke/lib/lib/*.* lib/lib

.FORCE:

.cc.o:
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

.c.o:
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
	@$(MKDIR) -p $(PREFIX)/lib
	@$(MKDIR) -p $(PREFIX)/include
	@$(MKDIR) -p $(PREFIX)/include/Clothoids
	@cp lib/lib/lib$(LIB_CLOTHOID).*              $(PREFIX)/lib
	@$(LDCONFIG) $(PREFIX)/lib
	@cp -f src/*.hh                               $(PREFIX)/include
	@cp -f src/Clothoids/*.hxx                    $(PREFIX)/include/Clothoids

install_all: install
	@cp lib/lib/lib$(LIB_CLOTHOID).*              $(PREFIX)/lib
	@cp submodules/Utils/lib/lib/*.*              $(PREFIX)/lib
	@cp submodules/quarticRootsFlocke/lib/lib/*.* $(PREFIX)/lib
	@cp -f submodules/quarticRootsFlocke/src/*.hh $(PREFIX)/include
	@cp -f submodules/Utils/src/*.hh              $(PREFIX)/include
	@cp -rf submodules/Utils/src/Utils            $(PREFIX)/include

install_as_framework: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	@cp -f src/*.hh                               $(PREFIX)/include/$(FRAMEWORK)
	@cp -f src/Clothoids/*.hxx                    $(PREFIX)/include/$(FRAMEWORK)/Clothoids
	@cp lib/lib/lib$(LIB_CLOTHOID)                $(PREFIX)/lib

install_as_framework_all: install_as_framework
	@cp -f submodules/quarticRootsFlocke/src/*.hh $(PREFIX)/include/$(FRAMEWORK)
	@cp -f submodules/Utils/src/*.hh              $(PREFIX)/include/$(FRAMEWORK)
	@cp -rf submodules/Utils/src/Utils            $(PREFIX)/include/$(FRAMEWORK)
	@cp lib/lib/lib$(LIB_CLOTHOID)                $(PREFIX)/lib
	@cp submodules/Utils/lib/lib/*.*              $(PREFIX)/lib
	@cp submodules/quarticRootsFlocke/lib/lib/*.* $(PREFIX)/lib

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
	rm -rf bin lib/*

depend:
	makedepend -- $(INC) $(CXXFLAGS) $(DEFS) -- $(SRCS)
# DO NOT DELETE

src/AABBtree.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/AABBtree.o: submodules/Utils/src/Utils/Utils.hxx
src/AABBtree.o: submodules/Utils/src/Utils/fmt/printf.h
src/AABBtree.o: submodules/Utils/src/Utils/fmt/ostream.h
src/AABBtree.o: submodules/Utils/src/Utils/fmt/format.h
src/AABBtree.o: submodules/Utils/src/Utils/fmt/core.h
src/AABBtree.o: submodules/Utils/src/Utils/fmt/chrono.h
src/AABBtree.o: submodules/Utils/src/Utils/fmt/locale.h
src/AABBtree.o: submodules/Utils/src/Utils/fmt/ostream.h
src/AABBtree.o: submodules/Utils/src/Utils/rang.hxx
src/AABBtree.o: submodules/Utils/src/Utils/Trace.hxx
src/AABBtree.o: submodules/Utils/src/Utils/Console.hxx
src/AABBtree.o: submodules/Utils/src/Utils/Malloc.hxx
src/AABBtree.o: submodules/Utils/src/Utils/Numbers.hxx
src/AABBtree.o: submodules/Utils/src/Utils/TicToc.hxx
src/AABBtree.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/AABBtree.o: submodules/Utils/src/Utils/Quaternion.hxx
src/AABBtree.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/AABBtree.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/AABBtree.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/AABBtree.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/AABBtree.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/AABBtree.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/AABBtree.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/AABBtree.o: src/Clothoids/ClothoidAsyPlot.hxx
src/Biarc.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/Biarc.o: submodules/Utils/src/Utils/Utils.hxx
src/Biarc.o: submodules/Utils/src/Utils/fmt/printf.h
src/Biarc.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Biarc.o: submodules/Utils/src/Utils/fmt/format.h
src/Biarc.o: submodules/Utils/src/Utils/fmt/core.h
src/Biarc.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Biarc.o: submodules/Utils/src/Utils/fmt/locale.h
src/Biarc.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Biarc.o: submodules/Utils/src/Utils/rang.hxx
src/Biarc.o: submodules/Utils/src/Utils/Trace.hxx
src/Biarc.o: submodules/Utils/src/Utils/Console.hxx
src/Biarc.o: submodules/Utils/src/Utils/Malloc.hxx
src/Biarc.o: submodules/Utils/src/Utils/Numbers.hxx
src/Biarc.o: submodules/Utils/src/Utils/TicToc.hxx
src/Biarc.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Biarc.o: submodules/Utils/src/Utils/Quaternion.hxx
src/Biarc.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/Biarc.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/Biarc.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/Biarc.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/Biarc.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/Biarc.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/Biarc.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/Biarc.o: src/Clothoids/ClothoidAsyPlot.hxx
src/BiarcList.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/BiarcList.o: submodules/Utils/src/Utils/Utils.hxx
src/BiarcList.o: submodules/Utils/src/Utils/fmt/printf.h
src/BiarcList.o: submodules/Utils/src/Utils/fmt/ostream.h
src/BiarcList.o: submodules/Utils/src/Utils/fmt/format.h
src/BiarcList.o: submodules/Utils/src/Utils/fmt/core.h
src/BiarcList.o: submodules/Utils/src/Utils/fmt/chrono.h
src/BiarcList.o: submodules/Utils/src/Utils/fmt/locale.h
src/BiarcList.o: submodules/Utils/src/Utils/fmt/ostream.h
src/BiarcList.o: submodules/Utils/src/Utils/rang.hxx
src/BiarcList.o: submodules/Utils/src/Utils/Trace.hxx
src/BiarcList.o: submodules/Utils/src/Utils/Console.hxx
src/BiarcList.o: submodules/Utils/src/Utils/Malloc.hxx
src/BiarcList.o: submodules/Utils/src/Utils/Numbers.hxx
src/BiarcList.o: submodules/Utils/src/Utils/TicToc.hxx
src/BiarcList.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/BiarcList.o: submodules/Utils/src/Utils/Quaternion.hxx
src/BiarcList.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/BiarcList.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/BiarcList.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/BiarcList.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/BiarcList.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/BiarcList.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/BiarcList.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/BiarcList.o: src/Clothoids/ClothoidAsyPlot.hxx
src/Circle.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/Circle.o: submodules/Utils/src/Utils/Utils.hxx
src/Circle.o: submodules/Utils/src/Utils/fmt/printf.h
src/Circle.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Circle.o: submodules/Utils/src/Utils/fmt/format.h
src/Circle.o: submodules/Utils/src/Utils/fmt/core.h
src/Circle.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Circle.o: submodules/Utils/src/Utils/fmt/locale.h
src/Circle.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Circle.o: submodules/Utils/src/Utils/rang.hxx
src/Circle.o: submodules/Utils/src/Utils/Trace.hxx
src/Circle.o: submodules/Utils/src/Utils/Console.hxx
src/Circle.o: submodules/Utils/src/Utils/Malloc.hxx
src/Circle.o: submodules/Utils/src/Utils/Numbers.hxx
src/Circle.o: submodules/Utils/src/Utils/TicToc.hxx
src/Circle.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Circle.o: submodules/Utils/src/Utils/Quaternion.hxx
src/Circle.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/Circle.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/Circle.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/Circle.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/Circle.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/Circle.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/Circle.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/Circle.o: src/Clothoids/ClothoidAsyPlot.hxx
src/Clothoid.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/Clothoid.o: submodules/Utils/src/Utils/Utils.hxx
src/Clothoid.o: submodules/Utils/src/Utils/fmt/printf.h
src/Clothoid.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Clothoid.o: submodules/Utils/src/Utils/fmt/format.h
src/Clothoid.o: submodules/Utils/src/Utils/fmt/core.h
src/Clothoid.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Clothoid.o: submodules/Utils/src/Utils/fmt/locale.h
src/Clothoid.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Clothoid.o: submodules/Utils/src/Utils/rang.hxx
src/Clothoid.o: submodules/Utils/src/Utils/Trace.hxx
src/Clothoid.o: submodules/Utils/src/Utils/Console.hxx
src/Clothoid.o: submodules/Utils/src/Utils/Malloc.hxx
src/Clothoid.o: submodules/Utils/src/Utils/Numbers.hxx
src/Clothoid.o: submodules/Utils/src/Utils/TicToc.hxx
src/Clothoid.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Clothoid.o: submodules/Utils/src/Utils/Quaternion.hxx
src/Clothoid.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/Clothoid.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/Clothoid.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/Clothoid.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/Clothoid.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/Clothoid.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/Clothoid.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/Clothoid.o: src/Clothoids/ClothoidAsyPlot.hxx
src/ClothoidAsyPlot.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/Utils.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/fmt/printf.h
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/fmt/ostream.h
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/fmt/format.h
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/fmt/core.h
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/fmt/chrono.h
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/fmt/locale.h
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/fmt/ostream.h
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/rang.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/Trace.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/Console.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/Malloc.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/Numbers.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/TicToc.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/Quaternion.hxx
src/ClothoidAsyPlot.o: submodules/Utils/src/Utils/Table.hxx
src/ClothoidAsyPlot.o: src/Clothoids/G2lib.hxx src/Clothoids/Triangle2D.hxx
src/ClothoidAsyPlot.o: src/Clothoids/BaseCurve.hxx src/Clothoids/AABBtree.hxx
src/ClothoidAsyPlot.o: src/Clothoids/Fresnel.hxx src/Clothoids/Line.hxx
src/ClothoidAsyPlot.o: src/Clothoids/BaseCurve_using.hxx
src/ClothoidAsyPlot.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/ClothoidAsyPlot.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/ClothoidAsyPlot.o: src/Clothoids/BiarcList.hxx
src/ClothoidAsyPlot.o: src/Clothoids/ClothoidList.hxx
src/ClothoidAsyPlot.o: src/Clothoids/ClothoidAsyPlot.hxx
src/ClothoidDistance.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/ClothoidDistance.o: submodules/Utils/src/Utils/Utils.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/fmt/printf.h
src/ClothoidDistance.o: submodules/Utils/src/Utils/fmt/ostream.h
src/ClothoidDistance.o: submodules/Utils/src/Utils/fmt/format.h
src/ClothoidDistance.o: submodules/Utils/src/Utils/fmt/core.h
src/ClothoidDistance.o: submodules/Utils/src/Utils/fmt/chrono.h
src/ClothoidDistance.o: submodules/Utils/src/Utils/fmt/locale.h
src/ClothoidDistance.o: submodules/Utils/src/Utils/fmt/ostream.h
src/ClothoidDistance.o: submodules/Utils/src/Utils/rang.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/Trace.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/Console.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/Malloc.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/Numbers.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/TicToc.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/Quaternion.hxx
src/ClothoidDistance.o: submodules/Utils/src/Utils/Table.hxx
src/ClothoidDistance.o: src/Clothoids/G2lib.hxx src/Clothoids/Triangle2D.hxx
src/ClothoidDistance.o: src/Clothoids/BaseCurve.hxx
src/ClothoidDistance.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/ClothoidDistance.o: src/Clothoids/Line.hxx
src/ClothoidDistance.o: src/Clothoids/BaseCurve_using.hxx
src/ClothoidDistance.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/ClothoidDistance.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/ClothoidDistance.o: src/Clothoids/BiarcList.hxx
src/ClothoidDistance.o: src/Clothoids/ClothoidList.hxx
src/ClothoidDistance.o: src/Clothoids/ClothoidAsyPlot.hxx
src/ClothoidG2.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/ClothoidG2.o: submodules/Utils/src/Utils/Utils.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/fmt/printf.h
src/ClothoidG2.o: submodules/Utils/src/Utils/fmt/ostream.h
src/ClothoidG2.o: submodules/Utils/src/Utils/fmt/format.h
src/ClothoidG2.o: submodules/Utils/src/Utils/fmt/core.h
src/ClothoidG2.o: submodules/Utils/src/Utils/fmt/chrono.h
src/ClothoidG2.o: submodules/Utils/src/Utils/fmt/locale.h
src/ClothoidG2.o: submodules/Utils/src/Utils/fmt/ostream.h
src/ClothoidG2.o: submodules/Utils/src/Utils/rang.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/Trace.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/Console.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/Malloc.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/Numbers.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/TicToc.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/Quaternion.hxx
src/ClothoidG2.o: submodules/Utils/src/Utils/Table.hxx
src/ClothoidG2.o: src/Clothoids/G2lib.hxx src/Clothoids/Triangle2D.hxx
src/ClothoidG2.o: src/Clothoids/BaseCurve.hxx src/Clothoids/AABBtree.hxx
src/ClothoidG2.o: src/Clothoids/Fresnel.hxx src/Clothoids/Line.hxx
src/ClothoidG2.o: src/Clothoids/BaseCurve_using.hxx src/Clothoids/Circle.hxx
src/ClothoidG2.o: src/Clothoids/Biarc.hxx src/Clothoids/Clothoid.hxx
src/ClothoidG2.o: src/Clothoids/PolyLine.hxx src/Clothoids/BiarcList.hxx
src/ClothoidG2.o: src/Clothoids/ClothoidList.hxx
src/ClothoidG2.o: src/Clothoids/ClothoidAsyPlot.hxx
src/ClothoidList.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/ClothoidList.o: submodules/Utils/src/Utils/Utils.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/fmt/printf.h
src/ClothoidList.o: submodules/Utils/src/Utils/fmt/ostream.h
src/ClothoidList.o: submodules/Utils/src/Utils/fmt/format.h
src/ClothoidList.o: submodules/Utils/src/Utils/fmt/core.h
src/ClothoidList.o: submodules/Utils/src/Utils/fmt/chrono.h
src/ClothoidList.o: submodules/Utils/src/Utils/fmt/locale.h
src/ClothoidList.o: submodules/Utils/src/Utils/fmt/ostream.h
src/ClothoidList.o: submodules/Utils/src/Utils/rang.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/Trace.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/Console.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/Malloc.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/Numbers.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/TicToc.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/Quaternion.hxx
src/ClothoidList.o: submodules/Utils/src/Utils/Table.hxx
src/ClothoidList.o: src/Clothoids/G2lib.hxx src/Clothoids/Triangle2D.hxx
src/ClothoidList.o: src/Clothoids/BaseCurve.hxx src/Clothoids/AABBtree.hxx
src/ClothoidList.o: src/Clothoids/Fresnel.hxx src/Clothoids/Line.hxx
src/ClothoidList.o: src/Clothoids/BaseCurve_using.hxx
src/ClothoidList.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/ClothoidList.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/ClothoidList.o: src/Clothoids/BiarcList.hxx
src/ClothoidList.o: src/Clothoids/ClothoidList.hxx
src/ClothoidList.o: src/Clothoids/ClothoidAsyPlot.hxx
src/Fresnel.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/Fresnel.o: submodules/Utils/src/Utils/Utils.hxx
src/Fresnel.o: submodules/Utils/src/Utils/fmt/printf.h
src/Fresnel.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Fresnel.o: submodules/Utils/src/Utils/fmt/format.h
src/Fresnel.o: submodules/Utils/src/Utils/fmt/core.h
src/Fresnel.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Fresnel.o: submodules/Utils/src/Utils/fmt/locale.h
src/Fresnel.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Fresnel.o: submodules/Utils/src/Utils/rang.hxx
src/Fresnel.o: submodules/Utils/src/Utils/Trace.hxx
src/Fresnel.o: submodules/Utils/src/Utils/Console.hxx
src/Fresnel.o: submodules/Utils/src/Utils/Malloc.hxx
src/Fresnel.o: submodules/Utils/src/Utils/Numbers.hxx
src/Fresnel.o: submodules/Utils/src/Utils/TicToc.hxx
src/Fresnel.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Fresnel.o: submodules/Utils/src/Utils/Quaternion.hxx
src/Fresnel.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/Fresnel.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/Fresnel.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/Fresnel.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/Fresnel.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/Fresnel.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/Fresnel.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/Fresnel.o: src/Clothoids/ClothoidAsyPlot.hxx
src/Fresnel.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
src/G2lib.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/G2lib.o: submodules/Utils/src/Utils/Utils.hxx
src/G2lib.o: submodules/Utils/src/Utils/fmt/printf.h
src/G2lib.o: submodules/Utils/src/Utils/fmt/ostream.h
src/G2lib.o: submodules/Utils/src/Utils/fmt/format.h
src/G2lib.o: submodules/Utils/src/Utils/fmt/core.h
src/G2lib.o: submodules/Utils/src/Utils/fmt/chrono.h
src/G2lib.o: submodules/Utils/src/Utils/fmt/locale.h
src/G2lib.o: submodules/Utils/src/Utils/fmt/ostream.h
src/G2lib.o: submodules/Utils/src/Utils/rang.hxx
src/G2lib.o: submodules/Utils/src/Utils/Trace.hxx
src/G2lib.o: submodules/Utils/src/Utils/Console.hxx
src/G2lib.o: submodules/Utils/src/Utils/Malloc.hxx
src/G2lib.o: submodules/Utils/src/Utils/Numbers.hxx
src/G2lib.o: submodules/Utils/src/Utils/TicToc.hxx
src/G2lib.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/G2lib.o: submodules/Utils/src/Utils/Quaternion.hxx
src/G2lib.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/G2lib.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/G2lib.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/G2lib.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/G2lib.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/G2lib.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/G2lib.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/G2lib.o: src/Clothoids/ClothoidAsyPlot.hxx
src/G2lib.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
src/G2lib_intersect.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/G2lib_intersect.o: submodules/Utils/src/Utils/Utils.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/fmt/printf.h
src/G2lib_intersect.o: submodules/Utils/src/Utils/fmt/ostream.h
src/G2lib_intersect.o: submodules/Utils/src/Utils/fmt/format.h
src/G2lib_intersect.o: submodules/Utils/src/Utils/fmt/core.h
src/G2lib_intersect.o: submodules/Utils/src/Utils/fmt/chrono.h
src/G2lib_intersect.o: submodules/Utils/src/Utils/fmt/locale.h
src/G2lib_intersect.o: submodules/Utils/src/Utils/fmt/ostream.h
src/G2lib_intersect.o: submodules/Utils/src/Utils/rang.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/Trace.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/Console.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/Malloc.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/Numbers.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/TicToc.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/Quaternion.hxx
src/G2lib_intersect.o: submodules/Utils/src/Utils/Table.hxx
src/G2lib_intersect.o: src/Clothoids/G2lib.hxx src/Clothoids/Triangle2D.hxx
src/G2lib_intersect.o: src/Clothoids/BaseCurve.hxx src/Clothoids/AABBtree.hxx
src/G2lib_intersect.o: src/Clothoids/Fresnel.hxx src/Clothoids/Line.hxx
src/G2lib_intersect.o: src/Clothoids/BaseCurve_using.hxx
src/G2lib_intersect.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/G2lib_intersect.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/G2lib_intersect.o: src/Clothoids/BiarcList.hxx
src/G2lib_intersect.o: src/Clothoids/ClothoidList.hxx
src/G2lib_intersect.o: src/Clothoids/ClothoidAsyPlot.hxx
src/Line.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/Line.o: submodules/Utils/src/Utils/Utils.hxx
src/Line.o: submodules/Utils/src/Utils/fmt/printf.h
src/Line.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Line.o: submodules/Utils/src/Utils/fmt/format.h
src/Line.o: submodules/Utils/src/Utils/fmt/core.h
src/Line.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Line.o: submodules/Utils/src/Utils/fmt/locale.h
src/Line.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Line.o: submodules/Utils/src/Utils/rang.hxx
src/Line.o: submodules/Utils/src/Utils/Trace.hxx
src/Line.o: submodules/Utils/src/Utils/Console.hxx
src/Line.o: submodules/Utils/src/Utils/Malloc.hxx
src/Line.o: submodules/Utils/src/Utils/Numbers.hxx
src/Line.o: submodules/Utils/src/Utils/TicToc.hxx
src/Line.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Line.o: submodules/Utils/src/Utils/Quaternion.hxx
src/Line.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/Line.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/Line.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/Line.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/Line.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/Line.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/Line.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/Line.o: src/Clothoids/ClothoidAsyPlot.hxx
src/PolyLine.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/PolyLine.o: submodules/Utils/src/Utils/Utils.hxx
src/PolyLine.o: submodules/Utils/src/Utils/fmt/printf.h
src/PolyLine.o: submodules/Utils/src/Utils/fmt/ostream.h
src/PolyLine.o: submodules/Utils/src/Utils/fmt/format.h
src/PolyLine.o: submodules/Utils/src/Utils/fmt/core.h
src/PolyLine.o: submodules/Utils/src/Utils/fmt/chrono.h
src/PolyLine.o: submodules/Utils/src/Utils/fmt/locale.h
src/PolyLine.o: submodules/Utils/src/Utils/fmt/ostream.h
src/PolyLine.o: submodules/Utils/src/Utils/rang.hxx
src/PolyLine.o: submodules/Utils/src/Utils/Trace.hxx
src/PolyLine.o: submodules/Utils/src/Utils/Console.hxx
src/PolyLine.o: submodules/Utils/src/Utils/Malloc.hxx
src/PolyLine.o: submodules/Utils/src/Utils/Numbers.hxx
src/PolyLine.o: submodules/Utils/src/Utils/TicToc.hxx
src/PolyLine.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/PolyLine.o: submodules/Utils/src/Utils/Quaternion.hxx
src/PolyLine.o: submodules/Utils/src/Utils/Table.hxx src/Clothoids/G2lib.hxx
src/PolyLine.o: src/Clothoids/Triangle2D.hxx src/Clothoids/BaseCurve.hxx
src/PolyLine.o: src/Clothoids/AABBtree.hxx src/Clothoids/Fresnel.hxx
src/PolyLine.o: src/Clothoids/Line.hxx src/Clothoids/BaseCurve_using.hxx
src/PolyLine.o: src/Clothoids/Circle.hxx src/Clothoids/Biarc.hxx
src/PolyLine.o: src/Clothoids/Clothoid.hxx src/Clothoids/PolyLine.hxx
src/PolyLine.o: src/Clothoids/BiarcList.hxx src/Clothoids/ClothoidList.hxx
src/PolyLine.o: src/Clothoids/ClothoidAsyPlot.hxx
src/Triangle2D.o: src/Clothoids.hh submodules/Utils/src/Utils.hh
src/Triangle2D.o: submodules/Utils/src/Utils/Utils.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/fmt/printf.h
src/Triangle2D.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Triangle2D.o: submodules/Utils/src/Utils/fmt/format.h
src/Triangle2D.o: submodules/Utils/src/Utils/fmt/core.h
src/Triangle2D.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Triangle2D.o: submodules/Utils/src/Utils/fmt/locale.h
src/Triangle2D.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Triangle2D.o: submodules/Utils/src/Utils/rang.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/Trace.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/Console.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/Malloc.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/Numbers.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/TicToc.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/Quaternion.hxx
src/Triangle2D.o: submodules/Utils/src/Utils/Table.hxx
src/Triangle2D.o: src/Clothoids/G2lib.hxx src/Clothoids/Triangle2D.hxx
src/Triangle2D.o: src/Clothoids/BaseCurve.hxx src/Clothoids/AABBtree.hxx
src/Triangle2D.o: src/Clothoids/Fresnel.hxx src/Clothoids/Line.hxx
src/Triangle2D.o: src/Clothoids/BaseCurve_using.hxx src/Clothoids/Circle.hxx
src/Triangle2D.o: src/Clothoids/Biarc.hxx src/Clothoids/Clothoid.hxx
src/Triangle2D.o: src/Clothoids/PolyLine.hxx src/Clothoids/BiarcList.hxx
src/Triangle2D.o: src/Clothoids/ClothoidList.hxx
src/Triangle2D.o: src/Clothoids/ClothoidAsyPlot.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/Utils/Table.hxx
submodules/Utils/src/CPUinfo.o: submodules/Utils/src/CPUinfo.hh
submodules/Utils/src/Console.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Table.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Table.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Table.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/Table.o: submodules/Utils/src/Utils/Table.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/TicToc.o: submodules/Utils/src/Utils/Table.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Table.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Table.hxx
submodules/Utils/src/fmt.o: submodules/Utils/src/Utils/fmt/os.cc
submodules/Utils/src/fmt.o: submodules/Utils/src/Utils/fmt/os.h
submodules/Utils/src/fmt.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/fmt.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/fmt.o: submodules/Utils/src/Utils/fmt/format.cc
submodules/Utils/src/fmt.o: submodules/Utils/src/Utils/fmt/format-inl.h
submodules/Utils/src/rang.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/rang.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/Quaternion.hxx
submodules/Utils/src/rang.o: submodules/Utils/src/Utils/Table.hxx
submodules/quarticRootsFlocke/src/PolynomialRoots-1-Quadratic.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-2-Cubic.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-3-Quartic.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-Jenkins-Traub.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-Utils.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-Utils.o: submodules/quarticRootsFlocke/src/PolynomialRoots-Utils.hh
submodules/GenericContainer/src/GenericContainer.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainer.o: submodules/GenericContainer/src/GenericContainerConfig.hh
submodules/GenericContainer/src/GenericContainerCinterface.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainerCinterface.o: submodules/GenericContainer/src/GenericContainerConfig.hh
submodules/GenericContainer/src/GenericContainerCinterface.o: submodules/GenericContainer/src/GenericContainerCinterface.h
submodules/GenericContainer/src/GenericContainerSerialize.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainerSerialize.o: submodules/GenericContainer/src/GenericContainerConfig.hh
submodules/GenericContainer/src/GenericContainerSerialize.o: submodules/GenericContainer/src/GenericContainerLibs.hh
submodules/GenericContainer/src/GenericContainerSupport.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainerSupport.o: submodules/GenericContainer/src/GenericContainerConfig.hh
submodules/GenericContainer/src/GenericContainerTables.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainerTables.o: submodules/GenericContainer/src/GenericContainerConfig.hh
