#pragma once

#include <geometry.hh>

class C0Coons {
public:
  using RationalCurve = std::pair<Geometry::PointVector, Geometry::DoubleVector>;

  // Constructors & evaluation
  C0Coons(const std::vector<RationalCurve> &boundaries);
  C0Coons(const std::vector<Geometry::PointVector> &boundaries); // (non-rational) Bezier curves
  Geometry::Point3D eval(const Geometry::Point2D &uv) const;
  Geometry::TriMesh eval(size_t resolution) const;

  // Mesh generation utilities
  Geometry::Point2DVector domain() const;
  Geometry::Point2DVector parameters(size_t resolution) const;
  Geometry::TriMesh meshTopology(size_t resolution) const;
  bool onEdge(size_t resolution, size_t index) const;

  // Rational Bezier curve evaluation utilities
  static Geometry::Point3D evalRational(const RationalCurve &curve, double u);
  static Geometry::Vector3D evalRationalDerivative(const RationalCurve &curve, double u);

private:
  void initialize();
  Geometry::Point3D evalRibbon(size_t i, const Geometry::Point2D &sd) const;

  size_t n_;
  std::vector<RationalCurve> boundaries_, opposites_;
  Geometry::Point2DVector domain_;
  Geometry::PointVector corners_;
};
