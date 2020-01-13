#include "c0coons.hh"

#include <algorithm>
#include <fstream>

using namespace Geometry;

std::vector<PointVector> readGBP(std::string filename) {
  std::ifstream f(filename.c_str());
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t n, d;
  f >> n >> d;
  size_t l = (d + 1) / 2;
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;          // # of control points

  Point3D p;
  f >> p[0] >> p[1] >> p[2];    // central control point

  std::vector<PointVector> result;
  PointVector ribbon;
  for (size_t i = 1, side = 0, col = 0; i < cp; ++i, ++col) {
    if (col > d) {
      col = 0;
      result.push_back(ribbon);
      auto last = ribbon.back();
      ribbon.clear();
      ribbon.push_back(last);
      if (++side >= n)
        break;
    }
    f >> p[0] >> p[1] >> p[2];
    ribbon.push_back(p);
  }

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

  auto curves = readGBP(argv[1]);
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
