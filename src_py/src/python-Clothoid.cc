/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-Clothoid.hh"
#include "python-ClothoidSpline-Interpolation.hh"
#include <pybind11/stl.h>
#include <stdexcept>

#ifdef _WIN32
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#endif

namespace G2lib {
  namespace python {
    void wrap_ClothoidCurve(py::module & m) {
      py::class_<ClothoidCurve, BaseCurve>(m, "ClothoidCurve")
        .def(py::init())
        .def(py::init<BaseCurve const &>())
        .def(py::init<ClothoidCurve const &>())
        .def(py::init<LineSegment const &>())
        .def(py::init<CircleArc const &>())
        .def(py::init<real_type, real_type, real_type, real_type, real_type, real_type>())
        .def(py::init([](std::tuple<real_type, real_type> P0, real_type theta0, 
                         std::tuple<real_type, real_type> P1, real_type theta1){
          real_type _P0[2] = {std::get<0>(P0), std::get<1>(P0)};
          real_type _P1[2] = {std::get<0>(P1), std::get<1>(P1)};
          return ClothoidCurve(_P0, theta0, _P1, theta1);
        }))
        .def("copy", &ClothoidCurve::copy)
        .def("build", py::overload_cast<real_type, real_type, real_type, real_type, real_type, real_type>(&ClothoidCurve::build))
        .def("build", py::overload_cast<LineSegment const &>(&ClothoidCurve::build))
        .def("build", py::overload_cast<CircleArc const &>(&ClothoidCurve::build))
        .def("build_forward", &ClothoidCurve::build_forward, 
             py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("kappa0"),
             py::arg("x1"), py::arg("y1"), py::arg("tol") = 1e-12)
        .def("build_G1", &ClothoidCurve::build_G1, 
             py::arg("x0"), py::arg("y0"), py::arg("theta0"),
             py::arg("x1"), py::arg("y1"), py::arg("theta1"), py::arg("tol") = 1e-12)
        .def("build_G1_D", [](ClothoidCurve * self,
                              real_type x0, real_type y0, real_type theta0,
                              real_type x1, real_type y1, real_type theta1,
                              std::tuple<real_type, real_type> L_D, std::tuple<real_type, real_type> k_D, 
                              std::tuple<real_type, real_type> dk_D, real_type tol = 1e-12) {
          real_type _L_D[2] = {std::get<0>(L_D), std::get<1>(L_D)};
          real_type _k_D[2] = {std::get<0>(k_D), std::get<1>(k_D)};
          real_type _dk_D[2] = {std::get<0>(dk_D), std::get<1>(dk_D)};
          return self->build_G1_D(x0, y0, theta0, x1, y1, theta1, _L_D, _k_D, _dk_D, tol);
        }, py::arg("x0"), py::arg("y0"), py::arg("theta0"), 
           py::arg("x1"), py::arg("y1"), py::arg("theta1"),
           py::arg("L_D"), py::arg("k_D"), py::arg("dk_D"), py::arg("tol") = 1e-12)
        .def("Pinfinity", [](ClothoidCurve * self, bool plus = true) {
          real_type x, y;
          self->Pinfinity(x, y, plus);
          return std::make_tuple(x, y);
        }, py::arg("plus") = true)
        .def("dkappa", &ClothoidCurve::dkappa)
        .def("thetaTotalVariation", &ClothoidCurve::thetaTotalVariation)
        .def("thetaMinMax", [](ClothoidCurve * self) {
          real_type th_min, th_max;
          self->thetaMinMax(th_min, th_max);
          return std::make_tuple(th_min, th_max);
        })
        .def("deltaTheta", &ClothoidCurve::deltaTheta)
        .def("curvatureMinMax", [](ClothoidCurve * self) {
          real_type cv_min, cv_max;
          self->curvatureMinMax(cv_min, cv_max);
          return std::make_tuple(cv_min, cv_max);
        })
        .def("curvatureTotalVariation", &ClothoidCurve::curvatureTotalVariation)
        .def("integralCurvature2", &ClothoidCurve::integralCurvature2)
        .def("integralJerk2", &ClothoidCurve::integralJerk2)
        .def("integralSnap2", &ClothoidCurve::integralSnap2)
        .def("optimized_sample_ISO", [](ClothoidCurve * self, real_type offs, int_type npts, real_type max_angle) {
          std::vector<real_type> ret;
          self->optimized_sample_ISO(offs, npts, max_angle, ret);
          return ret;
        })
        .def("optimized_sample_SAE", [](ClothoidCurve * self, real_type offs, int_type npts, real_type max_angle) {
          std::vector<real_type> ret;
          self->optimized_sample_SAE(offs, npts, max_angle, ret);
          return ret;
        })
        .def("closestPointBySample", [](ClothoidCurve * self, real_type ds, int_type qx, real_type qy) {
          real_type x, y, s;
          real_type v = self->closestPointBySample(ds, qx, qy, x, y, s);
          return std::make_tuple(v, std::make_tuple(x, y), s);
        })
        .def("distanceBySample", [](ClothoidCurve * self, real_type ds, int_type qx, real_type qy) {
          real_type s;
          real_type v = self->distanceBySample(ds, qx, qy, s);
          return std::make_tuple(v, s);
        })
        .def("distanceBySample", py::overload_cast<real_type, real_type, real_type>(&ClothoidCurve::distanceBySample, py::const_))
        .def("bbTriangle", [](ClothoidCurve * self, int_type icurve = 0) {
          Triangle2D t;
          bool v = self->bbTriangle(t, icurve);
          return std::make_tuple(v, t);
        }, py::arg("icurve") = 0)
        .def("bbTriangle_ISO", [](ClothoidCurve * self, real_type offs, int_type icurve = 0) {
          Triangle2D t;
          bool v = self->bbTriangle_ISO(offs, t, icurve);
          return std::make_tuple(v, t);
        }, py::arg("offs"), py::arg("icurve") = 0)
        .def("bbTriangle_SAE", [](ClothoidCurve * self, real_type offs, int_type icurve = 0) {
          Triangle2D t;
          bool v = self->bbTriangle_SAE(offs, t, icurve);
          return std::make_tuple(v, t);
        }, py::arg("offs"), py::arg("icurve") = 0)
        .def("bbTriangles", [](ClothoidCurve * self, real_type max_angle = m_pi/6, 
                                    real_type max_size  = 1e100, int_type icurve = 0) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles(tvec, max_angle, max_size, icurve);
          return tvec;
        }, py::arg("max_angle") = m_pi/6, py::arg("max_size") = 1e100, py::arg("icurve") = 0)
        .def("bbTriangles_ISO", [](ClothoidCurve * self, real_type offs, real_type max_angle, 
                                    real_type max_size, int_type icurve) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_ISO(offs, tvec, max_angle, max_size, icurve);
          return tvec;
        })
        .def("bbTriangles_SAE", [](ClothoidCurve * self, real_type offs, real_type max_angle, 
                                    real_type max_size, int_type icurve) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_SAE(offs, tvec, max_angle, max_size, icurve);
          return tvec;
        })
        .def("changeCurvilinearOrigin", &ClothoidCurve::changeCurvilinearOrigin)
        .def("build_AABBtree_ISO", &ClothoidCurve::build_AABBtree_ISO)
        .def("approximate_collision_ISO", &ClothoidCurve::approximate_collision_ISO);
    }

    void wrap_ClothoidSplineG2(py::module & m) {
      py::class_<ClothoidSplineG2>(m, "ClothoidSplineG2")
        .def(py::init<>())
        .def("setP1", &ClothoidSplineG2::setP1)
        .def("setP2", &ClothoidSplineG2::setP2)
        .def("setP3", &ClothoidSplineG2::setP3)
        .def("setP4", &ClothoidSplineG2::setP4)
        .def("setP5", &ClothoidSplineG2::setP5)
        .def("setP6", &ClothoidSplineG2::setP6)
        .def("setP7", &ClothoidSplineG2::setP7)
        .def("setP8", &ClothoidSplineG2::setP8)
        .def("setP9", &ClothoidSplineG2::setP9)
        .def("build", [](ClothoidSplineG2 * self, std::vector<real_type> x, std::vector<real_type> y) {
          int_type n = static_cast<int_type>(std::min(x.size(), y.size()));
          self->build(&x.front(), &y.front(), n);
        })
        .def("numPnts", &ClothoidSplineG2::numPnts)
        .def("numTheta", &ClothoidSplineG2::numTheta)
        .def("numConstraints", &ClothoidSplineG2::numConstraints)
        .def("guess", [](ClothoidSplineG2 * self) {
          size_t n = self->numPnts();
          std::vector<real_type> theta_guess(n), theta_min(n), theta_max(n);
          self->guess(&theta_guess.front(), &theta_min.front(), &theta_max.front());
          return std::make_tuple(theta_guess, theta_min, theta_max);
        })
        .def("objective", [](ClothoidSplineG2 * self, std::vector<real_type> theta) {
          real_type f;
          bool result = self->objective(&theta.front(), f);
          return std::make_tuple(result, f);
        })
        .def("gradient", [](ClothoidSplineG2 * self, std::vector<real_type> theta) {
          std::vector<real_type> g(self->numPnts());
          bool result = self->gradient(&theta.front(), &g.front());
          return std::make_tuple(result, g);
        })
        .def("constraints", [](ClothoidSplineG2 * self, std::vector<real_type> theta) {
          std::vector<real_type> c(self->numConstraints());
          bool result = self->constraints(&theta.front(), &c.front());
          return std::make_tuple(result, c);
        })
        .def("jacobian_nnz", &ClothoidSplineG2::jacobian_nnz)
        .def("jacobian_pattern", [](ClothoidSplineG2 * self) {
          std::vector<int_type> ii(self->jacobian_nnz()), jj(self->jacobian_nnz());
          bool result = self->jacobian_pattern(&ii.front(), &jj.front());
          return std::make_tuple(result, ii, jj);
        })
        .def("jacobian_pattern_matlab", [](ClothoidSplineG2 * self) {
          std::vector<real_type> ii(self->jacobian_nnz()), jj(self->jacobian_nnz());
          bool result = self->jacobian_pattern_matlab(&ii.front(), &jj.front());
          return std::make_tuple(result, ii, jj);
        })
        .def("jacobian", [](ClothoidSplineG2 * self, std::vector<real_type> theta) {
          std::vector<real_type> vals(self->jacobian_nnz());
          bool result = self->jacobian(&theta.front(), &vals.front());
          return std::make_tuple(result, vals);
        })
        .def("__str__", [](ClothoidSplineG2 * self) {
          std::ostringstream str;
          self->info(str);
          return str.str();
        });

      m.def("buildP1", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys, real_type theta_0, real_type theta_1) -> std::tuple<bool, G2lib::ClothoidList> {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator builder(xs, ys);
        bool info = builder.buildP1(theta_0, theta_1, result);
        return std::make_tuple(info, result);
      })
      .def("buildP2", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> std::tuple<bool, G2lib::ClothoidList> {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator builder(xs, ys);
        bool info = builder.buildP2(result);
        return std::make_tuple(info, result);
      })
      .def("buildP4", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> std::tuple<bool, G2lib::ClothoidList> {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator builder(xs, ys);
        bool info = builder.buildP4(result);
        return std::make_tuple(info, result);
      })
      .def("buildP5", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> std::tuple<bool, G2lib::ClothoidList> {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator builder(xs, ys);
        bool info = builder.buildP5(result);
        return std::make_tuple(info, result);
      })
      .def("buildP6", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> std::tuple<bool, G2lib::ClothoidList> {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator builder(xs, ys);
        bool info = builder.buildP6(result);
        return std::make_tuple(info, result);
      })
      .def("buildP7", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> std::tuple<bool, G2lib::ClothoidList> {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator builder(xs, ys);
        bool info = builder.buildP7(result);
        return std::make_tuple(info, result);
      })
      .def("buildP8", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> std::tuple<bool, G2lib::ClothoidList> {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator builder(xs, ys);
        bool info = builder.buildP8(result);
        return std::make_tuple(info, result);
      })
      .def("buildP9", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> std::tuple<bool, G2lib::ClothoidList> {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator builder(xs, ys);
        bool info = builder.buildP9(result);
        return std::make_tuple(info, result);
      });
    }

    void wrap_ClothoidList(py::module & m) {
      py::class_<ClothoidList, BaseCurve>(m, "ClothoidList")
        .def(py::init<>())
        .def(py::init<ClothoidList const &>())
        .def(py::init<LineSegment const &>())
        .def(py::init<CircleArc const &>())
        .def(py::init<Biarc const &>())
        .def(py::init<BiarcList const &>())
        .def(py::init<ClothoidCurve const &>())
        .def(py::init<PolyLine const &>())
        .def(py::init<BaseCurve const &>())
        .def("get", &ClothoidList::get)
        .def("init", &ClothoidList::init)
        .def("reserve", &ClothoidList::reserve)
        .def("push_back", py::overload_cast<LineSegment const &>(&ClothoidList::push_back))
        .def("push_back", py::overload_cast<CircleArc const &>(&ClothoidList::push_back))
        .def("push_back", py::overload_cast<Biarc const &>(&ClothoidList::push_back))
        .def("push_back", py::overload_cast<BiarcList const &>(&ClothoidList::push_back))
        .def("push_back", py::overload_cast<ClothoidCurve const &>(&ClothoidList::push_back))
        .def("push_back", py::overload_cast<PolyLine const &>(&ClothoidList::push_back))
        .def("push_back", py::overload_cast<real_type, real_type, real_type>(&ClothoidList::push_back))
        .def("push_back", py::overload_cast<real_type, real_type, real_type, real_type, real_type, real_type>(&ClothoidList::push_back))
        .def("push_back_G1", py::overload_cast<real_type, real_type, real_type>(&ClothoidList::push_back_G1))
        .def("push_back_G1", py::overload_cast<real_type, real_type, real_type, real_type, real_type, real_type>(&ClothoidList::push_back_G1))
        .def("build_G1", [](ClothoidList * self, std::vector<real_type> x, std::vector<real_type> y) {
          int_type n = static_cast<int_type>(std::min(x.size(), y.size()));
          return self->build_G1(n, &x.front(), &y.front());
        })
        .def("build_G1", [](ClothoidList * self, std::vector<real_type> x, std::vector<real_type> y, std::vector<real_type> theta) {
          int_type n = static_cast<int_type>(std::min({x.size(), y.size(), theta.size()}));
          return self->build_G1(n, &x.front(), &y.front(), &theta.front());
        })
        .def("build", [](ClothoidList * self, real_type x0, real_type y0, real_type theta0, std::vector<real_type> s, std::vector<real_type> kappa) {
          int_type n = static_cast<int_type>(std::min(s.size(), kappa.size()));
          return self->build(x0, y0, theta0, n, &s.front(), &kappa.front());
        })
        .def("getAtS", &ClothoidList::getAtS)
        .def("numSegment", &ClothoidList::numSegment)
        .def("findAtS", &ClothoidList::findAtS)
        .def("segment_length", &ClothoidList::segment_length)
        .def("segment_length_ISO", &ClothoidList::segment_length_ISO)
        .def("segment_length_SAE", &ClothoidList::segment_length_SAE)
        .def("bbTriangles", [](ClothoidList * self, real_type max_angle = m_pi/6, real_type max_size  = 1e100) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles(tvec, max_angle, max_size);
          return tvec;
        }, py::arg("max_angle") = m_pi/6, py::arg("max_size") = 1e100)
        .def("bbTriangles_ISO", [](ClothoidList * self, real_type offs, real_type max_angle, real_type max_size) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_ISO(offs, tvec, max_angle, max_size);
          return tvec;
        })
        .def("bbTriangles_SAE", [](ClothoidList * self, real_type offs, real_type max_angle = m_pi/18, real_type max_size  = 1e100) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_SAE(offs, tvec, max_angle, max_size);
          return tvec;
        }, py::arg("offs"), py::arg("max_angle") = m_pi/18, py::arg("max_size") = 1e100)
        .def("build_AABBtree_ISO", &ClothoidList::build_AABBtree_ISO)
        .def("getSK", [](ClothoidList * self) {
          int_type n = self->numSegment();
          std::vector<real_type> s(n), k(n);
          self->getSK(&s.front(), &k.front());
          return std::make_tuple(s, k);
        })
        .def("getSTK", [](ClothoidList * self) {
          int_type n = self->numSegment();
          std::vector<real_type> s(n), t(n), k(n);
          self->getSTK(&s.front(), &t.front(), &k.front());
          return std::make_tuple(s, t, k);
        })
        .def("getXY", [](ClothoidList * self) {
          int_type n = self->numSegment();
          std::vector<real_type> x(n), y(n);
          self->getXY(&x.front(), &y.front());
          return std::make_tuple(x, y);
        })
        .def("getDeltaTheta", [](ClothoidList * self) {
          int_type n = self->numSegment();
          std::vector<real_type> deltaTheta(n);
          self->getDeltaTheta(&deltaTheta.front());
          return deltaTheta;
        })
        .def("getDeltaKappa", [](ClothoidList * self) {
          int_type n = self->numSegment();
          std::vector<real_type> deltaKappa(n);
          self->getDeltaKappa(&deltaKappa.front());
          return deltaKappa;
        })
        .def("findST1", [](ClothoidList * self, real_type x, real_type y) {
          real_type s, t;
          self->findST1(x, y, s, t);
          return std::make_tuple(s, t);
        })
        .def("findST1", [](ClothoidList * self, int_type ibegin, int_type iend, real_type x, real_type y) {
          real_type s, t;
          self->findST1(ibegin, iend, x, y, s, t);
          return std::make_tuple(s, t);
        });
    }
  }
}