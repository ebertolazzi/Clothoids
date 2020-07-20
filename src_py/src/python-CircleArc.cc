/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-CircleArc.hh"
#include <pybind11/stl.h>

namespace G2lib {
  namespace python {
    void wrap_CircleArc(py::module & m) {
      py::class_<CircleArc, BaseCurve>(m, "CircleArc")
        .def(py::init<BaseCurve const &>())
        .def(py::init<LineSegment const &>())
        .def(py::init<CircleArc const &>())
        .def(py::init<real_type, real_type, real_type, real_type, real_type>())
        .def("copy", &CircleArc::copy)
        .def("build_G1", &CircleArc::build_G1)
        .def("build_3P", &CircleArc::build_3P)
        //.def("bbTriangle", [](CircleArc * self) {
        //  real_type x0, y0, x1, y1, x2, y2;
        //  self->bbTriangle(x0, y0, x1, y1, x2, y2);
        //  return std::make_tuple(
        //    std::make_tuple(x0, y0), 
        //    std::make_tuple(x1, y1), 
        //    std::make_tuple(x2, y2));
        //})
        .def("bbTriangle", [](CircleArc * self, real_type ss0 = 0, 
                                    real_type ss1 = 0, int_type icurve = 0) {
          Triangle2D t;
          bool v = self->bbTriangle(t, ss0, ss1, icurve);
          return std::make_tuple(v, t);
        }, py::arg("ss0") = 0, py::arg("ss1") = 0, py::arg("icurve") = 0)
        .def("bbTriangle_ISO", [](CircleArc * self, real_type offs, real_type ss0 = 0, 
                                  real_type ss1 = 0, int_type icurve = 0) {
          Triangle2D t;
          bool v = self->bbTriangle_ISO(offs, t, ss0, ss1, icurve);
          return std::make_tuple(v, t);
        }, py::arg("offs"), py::arg("ss0") = 0, py::arg("ss1") = 0, py::arg("icurve") = 0)
        .def("bbTriangle_SAE", [](CircleArc * self, real_type offs, real_type ss0 = 0, 
                                  real_type ss1 = 0, int_type icurve = 0) {
          Triangle2D t;
          bool v = self->bbTriangle_SAE(offs, t, ss0, ss1, icurve);
          return std::make_tuple(v, t);
        }, py::arg("offs"), py::arg("ss0") = 0, py::arg("ss1") = 0, py::arg("icurve") = 0)
        .def("bbTriangles", [](CircleArc * self, real_type max_angle = m_pi/18, 
                                    real_type max_size  = 1e100, int_type icurve = 0) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles(tvec, max_angle, max_size, icurve);
          return tvec;
        }, py::arg("max_angle") = m_pi/18, py::arg("max_size") = 1e100, py::arg("icurve") = 0)
        .def("bbTriangles_ISO", [](CircleArc * self, real_type offs, real_type max_angle, 
                                    real_type max_size, int_type icurve) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_ISO(offs, tvec, max_angle, max_size, icurve);
          return tvec;
        })
        .def("bbTriangles_SAE", [](CircleArc * self, real_type offs, real_type max_angle, 
                                    real_type max_size, int_type icurve) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_SAE(offs, tvec, max_angle, max_size, icurve);
          return tvec;
        })
        .def("collision", &CircleArc::collision)
        .def("intersect", [](CircleArc * self, CircleArc const & C, bool swap_s_vals) {
          IntersectList ilist;
          self->intersect(C, ilist, swap_s_vals);
          return ilist;
        })
        .def("collision_ISO", &CircleArc::collision_ISO)
        .def("intersect_ISO", [](CircleArc * self, real_type offs, CircleArc const & C, real_type Coffs, bool swap_s_vals) {
          IntersectList ilist;
          self->intersect_ISO(offs, C, Coffs, ilist, swap_s_vals);
          return ilist;
        })
        .def("sinTheta0", &CircleArc::sinTheta0)
        .def("cosTheta0", &CircleArc::cosTheta0)
        .def("curvature", &CircleArc::curvature)
        .def("lenTolerance", &CircleArc::lenTolerance)
        .def("delta_theta", &CircleArc::delta_theta)
        .def("thetaTotalVariation", &CircleArc::thetaTotalVariation)
        .def("thetaMinMax", [](CircleArc * self) {
          real_type th_min, th_max;
          self->thetaMinMax(th_min, th_max);
          return std::make_tuple(th_min, th_max);
        })
        .def("changeCurvilinearOrigin", &CircleArc::changeCurvilinearOrigin)
        .def("center", [](CircleArc * self) {
          real_type x, y;
          self->center(x, y);
          return std::make_tuple(x, y);
        })
        .def("ray", &CircleArc::ray)
        .def("paramNURBS", [](CircleArc * self) {
          int_type n_pnts, n_knots;
          self->paramNURBS(n_knots, n_pnts);
          return std::make_tuple(n_knots, n_pnts);
        })
        .def("toNURBS", [](CircleArc * self){
          using Point = real_type[3];
          int_type n_pnts, n_knots;
          self->paramNURBS(n_knots, n_pnts);
          real_type * knots = new real_type[n_knots];
          Point * Poly = new Point[n_pnts];
          self->toNURBS(knots, Poly);

          std::vector<std::tuple<real_type, real_type, real_type>> Poly_vec;
          for (int_type i = 0; i < n_pnts; i++) {
            Poly_vec.push_back(std::make_tuple(Poly[i][0], Poly[i][1], Poly[i][2]));
          }
          std::vector<real_type> knots_vec;
          for (int_type i = 0; i < n_knots; i++) {
            knots_vec.push_back(knots[i]);
          }
          auto ret = std::make_tuple(knots_vec, Poly_vec);
          delete[] Poly;
          delete[] knots;
          return ret;
        });
    }
  }
}