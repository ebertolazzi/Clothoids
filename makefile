# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

LIB_CLOTHOID = libClothoid.a

CC   = gcc
CXX  = g++
INC  = -I./src -I./include
LIBS = -L./lib -lClothoid
DEFS =

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  LIBS     = -static -L./lib -lClothoid
  CXXFLAGS = -Wall -O3 -fPIC -Wno-sign-compare
  AR       = ar rcs
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  CC       = clang
  CXX      = clang++
  #CC       = gcc-7
  #CXX      = g++-7
  LIBS     = -L./lib -lClothoid
  CXXFLAGS = -Wall -O3 -fPIC -Wno-sign-compare
  AR       = libtool -static -o
endif

SRCS = \
src/Biarc.cc \
src/Circle.cc \
src/ClothoidAsyPlot.cc \
src/ClothoidG1.cc \
src/ClothoidG2.cc \
src/CubicRootsFlocke.cc \
src/Fresnel.cc \
src/G2lib.cc \
src/Line.cc \
src/Triangle2D.cc

OBJS  = $(SRCS:.cc=.o)
DEPS  = src/Clothoid.hh src/CubicRootsFlocke.hh
MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Clothoid

all: lib
	@$(MKDIR) bin
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2         tests-cpp/testG2.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2stat     tests-cpp/testG2stat.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2stat2arc tests-cpp/testG2stat2arc.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2statCLC  tests-cpp/testG2statCLC.cc $(LIBS)
	#$(CXX) $(INC) $(CXXFLAGS) -o bin/testG2plot tests-cpp/testG2plot.cc $(LIBS)

lib: lib/$(LIB_CLOTHOID)

src/%.o: src/%.cc $(DEPS)
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@ 

src/%.o: src/%.c $(DEPS)
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

lib/libClothoid.a: $(OBJS)
	@$(MKDIR) lib
	$(AR) lib/libClothoid.a $(OBJS) 

lib/libClothoid.dylib: $(OBJS)
	@$(MKDIR) lib
	$(CXX) -shared -o lib/libClothoid.dylib $(OBJS) 

lib/libClothoid.so: $(OBJS)
	@$(MKDIR) lib
	$(CXX) -shared -o lib/libClothoid.so $(OBJS) 

install: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include
	cp src/Clothoid.hh         $(PREFIX)/include
	cp src/CubicRootsFlocke.hh $(PREFIX)/include
	cp lib/$(LIB_CLOTHOID)     $(PREFIX)/lib

install_as_framework: lib
	@$(MKDIR) $(PREFIX)/lib
	@$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp src/Clothoid.hh         $(PREFIX)/include/$(FRAMEWORK)
	cp src/CubicRootsFlocke.hh $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_CLOTHOID)     $(PREFIX)/lib

run:
	./bin/testG2
	./bin/testG2stat
	./bin/testG2stat2arc
	./bin/testG2statCLC

doc:
	doxygen
	
clean:
	rm -f lib/libClothoid.* lib/libClothoid.* src/*.o

	rm -rf bin
	