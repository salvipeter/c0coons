all: example

INCLUDES=-I../../libgeom
LIBS=-L../../libgeom/release -lgeom
CXXFLAGS=-Wall -pedantic -std=c++17 $(INCLUDES)

example: example.o c0coons.o curves.o
	g++ -o $@ $^ $(LIBS)

c0coons.o: c0coons.cc c0coons.hh

curves.o: curves.cc curves.hh
