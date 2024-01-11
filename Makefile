all: libc0coons.a example

INCLUDES=-I../libgeom
LIBS=-L../libgeom/release -lgeom -L. -lc0coons
CXXFLAGS=-Wall -pedantic -std=c++17 $(INCLUDES)

example: example.o libc0coons.a
	g++ -o $@ $< $(LIBS)

libc0coons.a: c0coons.o curves.o
	$(AR) r -o $@ $^

c0coons.o: c0coons.cc c0coons.hh

curves.o: curves.cc curves.hh
