#include "c0coons.hh"

#include <algorithm>
#include <fstream>

using namespace Geometry;

std::vector<PointVector> readBezierBoundary(std::string filename) {
  std::ifstream f(filename.c_str());
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t n, d;
  f >> n >> d;

  Point3D p;
  f >> p[0] >> p[1] >> p[2];    // central control point (not used)

  std::vector<PointVector> result;
  PointVector ribbon;
  for (size_t col = 0, side = 0; true; ++col) {
    if (col > d) {
      col = 1;
      result.push_back(ribbon);
      auto last = ribbon.back();
      ribbon.clear();
      ribbon.push_back(last);
      if (++side == n)
        break;
    }
    f >> p[0] >> p[1] >> p[2];
    ribbon.push_back(p);
  }
  result.back().back() = result.front().front();

  return result;
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <model.gbp> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 15;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  C0Coons surface(readBezierBoundary(argv[1]));
  surface.eval(resolution).writeOBJ("test.obj");
}
