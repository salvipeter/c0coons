#pragma once

#include <geometry.hh>

namespace Geometry {

  class CurveType {             // base type, parameterized in [0,1]
  public:
    virtual ~CurveType() { }
    virtual Point3D eval(double u) const = 0;
    virtual Vector3D evalDerivative(double u) const = 0;
  };

  class RationalBezierCurve : public CurveType {
  public:
    RationalBezierCurve(const PointVector &cpts, const DoubleVector &weights);
    virtual Point3D eval(double u) const override;
    virtual Vector3D evalDerivative(double u) const override;
  private:
    PointVector cp;
    DoubleVector wi;
  };

  class BezierCurve : public RationalBezierCurve {
  public:
    BezierCurve(const PointVector &cpts);
  };

  class OnePointCurve : public CurveType { // degenerate curve having only 1 point
  public:
    OnePointCurve(const Point3D &p) : p(p) { }
    virtual Point3D eval(double) const override { return p; }
    virtual Vector3D evalDerivative(double) const override { return { 0, 0, 0 }; }
  private:
    Point3D p;
  };

  class BSplineCurve : public CurveType { // just a wrapper on BSCurve
  public:
    BSplineCurve(const BSCurve &curve) : c(curve) { c.normalize(); }
    virtual Point3D eval(double u) const override { return c.eval(u); }
    virtual Vector3D evalDerivative(double u) const override {
      VectorVector der;
      c.eval(u, 1, der);
      return der[1];
    }
  private:
    BSCurve c;
  };

} // namespace Geometry
