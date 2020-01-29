#include "curves.hh"

namespace Geometry {

static DoubleVector bernstein(size_t n, double u) {
  DoubleVector coeff; coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double  tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
  return coeff;
}

RationalBezierCurve::RationalBezierCurve(const PointVector &cpts, const DoubleVector &weights)
  : cp(cpts), wi(weights)
{
}

Point3D
RationalBezierCurve::eval(double u) {
  size_t n = cp.size() - 1;
  DoubleVector coeff = bernstein(n, u);
  Point3D p(0.0, 0.0, 0.0);
  double w = 0;
  for (size_t k = 0; k <= n; ++k) {
    p += cp[k] * wi[k] * coeff[k];
    w += wi[k] * coeff[k];
  }
  return p / w;
}

Vector3D
RationalBezierCurve::evalDerivative(double u) {
  size_t n = cp.size() - 1;

  PointVector dcp;
  DoubleVector dwi;
  for (size_t i = 0; i < n; ++i) {
    dcp.push_back((cp[i+1] * wi[i+1] - cp[i] * wi[i]) * n);
    dwi.push_back((wi[i+1] - wi[i]) * n);
  }

  DoubleVector coeff = bernstein(n, u);
  Point3D p(0.0, 0.0, 0.0);
  double w = 0;
  for (size_t k = 0; k <= n; ++k) {
    p += cp[k] * wi[k] * coeff[k];
    w += wi[k] * coeff[k];
  }

  coeff = bernstein(n - 1, u);
  Point3D dp(0.0, 0.0, 0.0);
  double dw = 0;
  for (size_t k = 0; k < n; ++k) {
    dp += dcp[k] * coeff[k];
    dw += dwi[k] * coeff[k];
  }

  return (dp - p / w * dw) / w;
}

BezierCurve::BezierCurve(const PointVector &cpts)
  : RationalBezierCurve(cpts, DoubleVector(cpts.size(), 1))
{
}

} // namespace Geometry
