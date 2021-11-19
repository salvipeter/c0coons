project(c0coons)

set(LIBGEOM_DIR "../geom" CACHE FILEPATH "Libgeom path" FORCE)

include_directories("${LIBGEOM_DIR}")

add_library(c0coons STATIC
	c0coons.cc
	c0coons.hh
	curves.cc
	curves.hh
)
