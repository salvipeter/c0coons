#include "c0coons.hh"

#include <algorithm>
#include <fstream>

using namespace Geometry;

std::vector<std::shared_ptr<CurveType>> readRationalBoundaries(std::string filename) {
  std::ifstream f(filename.c_str());
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t n, d;
  f >> n >> d;

  std::vector<std::shared_ptr<CurveType>> result;
  PointVector points(d + 1);
  DoubleVector weights(d + 1);

  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= d; ++j)
      f >> points[j][0] >> points[j][1] >> points[j][2] >> weights[j];
    auto curve = std::make_shared<RationalBezierCurve>(points, weights);
    result.push_back(std::dynamic_pointer_cast<CurveType>(curve));
  }

  return result;
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <model.rbl> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 15;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  C0Coons surface(readRationalBoundaries(argv[1]));
  surface.eval(resolution).writeOBJ("test.obj");
}
