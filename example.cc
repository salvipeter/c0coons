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

  auto curves = readBezierBoundary(argv[1]);
  std::vector<C0Coons::RationalCurve> boundaries;
  std::transform(curves.begin(), curves.end(), std::back_inserter(boundaries),
                 [&](const PointVector &points) {
                   return std::make_pair(points, DoubleVector(points.size(), 1));
                 });

  C0Coons surface(boundaries);

  TriMesh mesh = surface.meshTopology(resolution);
  Point2DVector uvs = surface.parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) { return surface.eval(uv); });
  mesh.setPoints(points);
  mesh.writeOBJ("test.obj");
}
