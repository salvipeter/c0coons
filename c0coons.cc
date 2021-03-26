#include "c0coons.hh"

#include <algorithm>
#include <cmath>
#include <numeric>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Geometry;

static double inrange(double min, double x, double max) {
  if (x < min)
    return min;
  if (x > max)
    return max;
  return x;
}

void
C0Coons::initialize() {         // assumes that n_ and boundaries_ are already set
  // Setup Domain
  if (n_ == 4)
    domain_ = { {1,1}, {-1,1}, {-1,-1}, {1,-1} };
  else {
    double alpha = 2.0 * M_PI / n_;
    for (size_t i = 0; i < n_; ++i)
      domain_.emplace_back(std::cos(alpha * i), std::sin(alpha * i));
  }

  // Store corners
  for (const auto &b : boundaries_)
    corners_.push_back(b->eval(0));

  // Setup opposite curves
  if (n_ == 3) { // 0-degree rational Bezier curves
    opposites_.push_back(std::dynamic_pointer_cast<CurveType>
                         (std::make_shared<OnePointCurve>(corners_[2])));
    opposites_.push_back(std::dynamic_pointer_cast<CurveType>
                         (std::make_shared<OnePointCurve>(corners_[0])));
    opposites_.push_back(std::dynamic_pointer_cast<CurveType>
                         (std::make_shared<OnePointCurve>(corners_[1])));
  } else
    for (size_t i = 0; i < n_; ++i) {
      size_t imm = (i + n_ - 2) % n_, im = (i + n_ - 1) % n_, ipp = (i + 2) % n_;
      auto p1 = corners_[im];
      auto p2 = p1 - boundaries_[imm]->evalDerivative(1) / 3;
      auto q1 = corners_[ipp];
      auto q2 = q1 + boundaries_[ipp]->evalDerivative(0) / 3;
      PointVector cpts = { q1, q2, p2, p1 };
      opposites_.push_back(std::dynamic_pointer_cast<CurveType>
                           (std::make_shared<BezierCurve>(cpts)));
    }
}

C0Coons::C0Coons(std::vector<std::shared_ptr<CurveType>> boundaries)
  : n_(boundaries.size()), boundaries_(std::move(boundaries))
{
  initialize();
}

Point3D
C0Coons::evalRibbon(size_t i, const Point2D &sd) const {
  size_t im = (i + n_ - 1) % n_, ip = (i + 1) % n_, ipp = (i + 2) % n_;
  auto s = inrange(0, sd[0], 1), d = inrange(0, sd[1], 1);
  auto s1 = inrange(0, 1 - s, 1), d1 = inrange(0, 1 - d, 1);
  auto p1 = boundaries_[i]->eval(s) * d1 + opposites_[i]->eval(s1) * d;
  auto p2 = boundaries_[im]->eval(d1) * s1 + boundaries_[ip]->eval(d) * s;
  auto p12 = (corners_[i] * s1 + corners_[ip] * s) * d1 +
    (corners_[im] * s1 + corners_[ipp] * s) * d;
  return p1 + p2 - p12;
}

static DoubleVector wachspress(const Point2DVector &domain, const Point2D &uv) {
  size_t n = domain.size();
  Vector2DVector vectors; vectors.reserve(n);
  std::transform(domain.begin(), domain.end(), std::back_inserter(vectors),
                 [uv](const Point2D &p) { return uv - p; });

  DoubleVector areas; areas.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    const Vector2D &si = vectors[i];
    const Vector2D &si1 = vectors[(i+1)%n];
    areas.push_back((si[0] * si1[1] - si[1] * si1[0]) / 2.0);
  }

  DoubleVector l; l.reserve(n);

  for (size_t i = 0; i < n; ++i) {
    size_t i_1 = (i + n - 1) % n, i1 = (i + 1) % n;
    double Ai = 1.0, Ai_1 = 1.0, Ai_1i = 1.0;
    for (size_t j = 0; j < n; ++j) {
      if (j == i)
        Ai_1 *= areas[j];
      else if (j == i_1)
        Ai *= areas[j];
      else {
        Ai_1 *= areas[j];
        Ai *= areas[j];
        Ai_1i *= areas[j];
      }
    }
    const Vector2D &si_1 = vectors[i_1];
    const Vector2D &si1 = vectors[i1];
    double Bi = (si_1[0] * si1[1] - si_1[1] * si1[0]) / 2.0;
    l.push_back(Ai_1 + Ai - Bi * Ai_1i);
  }

  double sum = std::accumulate(l.begin(), l.end(), 0.0);
  std::transform(l.begin(), l.end(), l.begin(), [sum](double x) { return x / sum; });
  return l;
}

Point3D
C0Coons::eval(const Point2D &uv) const {
  auto bc = wachspress(domain_, uv);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i) {
    double h1 = bc[(i+n_-1)%n_] + bc[i];
    double s = h1 > epsilon ? bc[i] / h1 : 0.5; // kutykurutty
    p += evalRibbon(i, { s, 1 - h1 }) * h1 / 2;
  }
  return p;
}

TriMesh
C0Coons::eval(size_t resolution) const {
  TriMesh mesh = meshTopology(resolution);
  Point2DVector uvs = parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) { return eval(uv); });
  mesh.setPoints(points);
  return mesh;
}

Point2DVector
C0Coons::domain() const {
  return domain_;
}

static size_t meshSize(size_t n, size_t resolution) {
  if (n == 3)
    return (resolution + 1) * (resolution + 2) / 2;
  if (n == 4)
    return (resolution + 1) * (resolution + 1);
  return 1 + n * resolution * (resolution + 1) / 2;
}

Point2DVector
C0Coons::parameters(size_t resolution) const {
  size_t size = meshSize(n_, resolution);
  Point2DVector result;
  result.reserve(size);

  if (n_ == 3) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = domain_[0] * u + domain_[2] * (1 - u);
      auto q = domain_[1] * u + domain_[2] * (1 - u);
      for (size_t k = 0; k <= j; ++k) {
        double v = j == 0 ? 1.0 : (double)k / j;
        result.push_back(p * (1 - v) + q * v);
      }
    }
  } else if (n_ == 4) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = domain_[0] * (1 - u) + domain_[1] * u;
      auto q = domain_[3] * (1 - u) + domain_[2] * u;
      for (size_t k = 0; k <= resolution; ++k) {
        double v = (double)k / resolution;
        result.push_back(p * (1 - v) + q * v);
      }
    }
  } else { // n_ > 4
    Point2D center(0.0, 0.0);
    result.push_back(center);
    for (size_t j = 1; j <= resolution; ++j) {
      double u = (double)j / (double)resolution;
      for (size_t k = 0; k < n_; ++k)
        for (size_t i = 0; i < j; ++i) {
          double v = (double)i / (double)j;
          Point2D ep = domain_[(k+n_-1)%n_] * (1.0 - v) + domain_[k] * v;
          Point2D p = center * (1.0 - u) + ep * u;
          result.push_back(p);
        }
    }
  }
  return result;
}

TriMesh
C0Coons::meshTopology(size_t resolution) const {
  TriMesh mesh;
  mesh.resizePoints(meshSize(n_, resolution));

  if (n_ == 3) {
    size_t prev = 0, current = 1;
    for (size_t i = 0; i < resolution; ++i) {
      for (size_t j = 0; j < i; ++j) {
        mesh.addTriangle(current + j, current + j + 1, prev + j);
        mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
      }
      mesh.addTriangle(current + i, current + i + 1, prev + i);
      prev = current;
      current += i + 2;
    }
  } else if (n_ == 4) {
    for (size_t i = 0; i < resolution; ++i)
      for (size_t j = 0; j < resolution; ++j) {
        size_t index = i * (resolution + 1) + j;
        mesh.addTriangle(index, index + resolution + 1, index + 1);
        mesh.addTriangle(index + 1, index + resolution + 1, index + resolution + 2);
      }
  } else { // n_ > 4
    size_t inner_start = 0, outer_vert = 1;
    for (size_t layer = 1; layer <= resolution; ++layer) {
      size_t inner_vert = inner_start, outer_start = outer_vert;
      for (size_t side = 0; side < n_; ++side) {
        size_t vert = 0;
        while(true) {
          size_t next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
          mesh.addTriangle(inner_vert, outer_vert, next_vert);
          ++outer_vert;
          if (++vert == layer)
            break;
          size_t inner_next = (side == n_ - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
          mesh.addTriangle(inner_vert, next_vert, inner_next);
          inner_vert = inner_next;
        }
      }
      inner_start = outer_start;
    }
  }
  return mesh;
}

bool
C0Coons::onEdge(size_t resolution, size_t index) const {
  if (n_ == 3) {
    if (index >= meshSize(3, resolution) - resolution - 1)
      return true;
    auto issquare = [](size_t n) {
                      size_t root = std::round(std::sqrt(n));
                      return root * root == n;
                    };
    size_t n = index * 8 + 1;
    return issquare(n) || issquare(n + 8);
  }
  if (n_ == 4) {
    return index <= resolution || index >= (resolution + 1) * resolution ||
      index % (resolution + 1) == 0 || index % (resolution + 1) == resolution;
  }
  return index >= meshSize(n_, resolution) - n_ * resolution;
}
