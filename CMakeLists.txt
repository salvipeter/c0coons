project(c0coons)

include_directories("../geom")

add_library(c0coons STATIC
	c0coons.cc
	c0coons.hh
	curves.cc
	curves.hh
)
