#pragma once

#include "curves.hh"

class C0Coons {
public:
  // Constructor
  explicit C0Coons(std::vector<std::shared_ptr<Geometry::CurveType>> boundaries);

  // Evaluation
  Geometry::Point3D eval(const Geometry::Point2D &uv) const;
  Geometry::TriMesh eval(size_t resolution) const;

  // Mesh generation utilities
  Geometry::Point2DVector domain() const;
  Geometry::Point2DVector parameters(size_t resolution) const;
  Geometry::TriMesh meshTopology(size_t resolution) const;
  bool onEdge(size_t resolution, size_t index) const;

private:
  void initialize();
  Geometry::Point3D evalRibbon(size_t i, const Geometry::Point2D &sd) const;

  size_t n_;
  std::vector<std::shared_ptr<Geometry::CurveType>> boundaries_, opposites_;
  Geometry::Point2DVector domain_;
  Geometry::PointVector corners_;
};
