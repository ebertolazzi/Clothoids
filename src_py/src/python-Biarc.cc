/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-Biarc.hh"
#include "pybind11/stl.h"

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
    void wrap_Biarc(py::module & m) {
      py::class_<Biarc, BaseCurve>(m, "Biarc")
        .def(py::init())
        .def(py::init<Biarc const &>())
        .def(py::init<real_type, real_type, real_type, real_type, real_type, real_type>())
        .def(py::init<BaseCurve const &>())
        .def("copy", &Biarc::copy)
        .def("build", &Biarc::build)
        .def("build_3P", &Biarc::build_3P)
        .def("getC0", &Biarc::C0)
        .def("getC1", &Biarc::C1)
        .def("xMiddle", &Biarc::xMiddle)
        .def("yMiddle", &Biarc::yMiddle)
        .def("thetaMiddle", &Biarc::thetaMiddle)
        .def("kappa0", &Biarc::kappa0)
        .def("kappa1", &Biarc::kappa1)
        .def("length0", &Biarc::length0)
        .def("length1", &Biarc::length1)
        .def("delta_theta", &Biarc::delta_theta)
        .def("bbTriangles", [](Biarc * self, real_type max_angle = Utils::m_pi/18, 
                                    real_type max_size  = 1e100, int_type icurve = 0) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles(tvec, max_angle, max_size, icurve);
          return tvec;
        }, py::arg("max_angle") = Utils::m_pi/18, py::arg("max_size") = 1e100, py::arg("icurve") = 0)
        .def("bbTriangles_ISO", [](Biarc * self, real_type offs, real_type max_angle = Utils::m_pi/18, 
                                    real_type max_size  = 1e100, int_type icurve = 0) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_ISO(offs, tvec, max_angle, max_size, icurve);
          return tvec;
        }, py::arg("offs"), py::arg("max_angle") = Utils::m_pi/18, py::arg("max_size") = 1e100, py::arg("icurve") = 0)
        .def("bbTriangles_SAE", [](Biarc * self, real_type offs, real_type max_angle = Utils::m_pi/18, 
                                    real_type max_size  = 1e100, int_type icurve = 0) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_SAE(offs, tvec, max_angle, max_size, icurve);
          return tvec;
        }, py::arg("offs"), py::arg("max_angle") = Utils::m_pi/18, py::arg("max_size") = 1e100, py::arg("icurve") = 0)
        .def("collision", &Biarc::collision)
        .def("collision_ISO", &Biarc::collision_ISO)
        .def("intersect", [](Biarc * self, Biarc const & B, bool swap_s_vals) {
          IntersectList ilist;
          self->intersect(B, ilist, swap_s_vals);
          return ilist;
        })
        .def("intersect_ISO", [](Biarc * self, real_type offs, Biarc const & B, real_type Coffs, bool swap_s_vals) {
          IntersectList ilist;
          self->intersect_ISO(offs, B, Coffs, ilist, swap_s_vals);
          return ilist;
        });
    }

    void wrap_BiarcList(py::module & m) {
      py::class_<BiarcList, BaseCurve>(m, "BiarcList")
        .def(py::init<>())
        .def(py::init<BaseCurve const &>())
        .def(py::init<BiarcList const &>())
        .def(py::init<PolyLine const &>())
        .def(py::init<LineSegment const &>())
        .def(py::init<CircleArc const &>())
        .def(py::init<Biarc const &>())
        .def("get", &BiarcList::get)
        .def("copy", &BiarcList::copy)
        .def("init", &BiarcList::init)
        .def("reserve", &BiarcList::reserve)
        .def("push_back", py::overload_cast<LineSegment const &>(&BiarcList::push_back))
        .def("push_back", py::overload_cast<CircleArc const &>(&BiarcList::push_back))
        .def("push_back", py::overload_cast<Biarc const &>(&BiarcList::push_back))
        .def("push_back", py::overload_cast<PolyLine const &>(&BiarcList::push_back))
        .def("push_back_G1", py::overload_cast<real_type, real_type, real_type>(&BiarcList::push_back_G1))
        .def("push_back_G1", py::overload_cast<real_type, real_type, real_type, real_type, real_type, real_type>(&BiarcList::push_back_G1))
        .def("build_G1", [](BiarcList * self, std::vector<real_type> x, std::vector<real_type> y) {
          int_type n = static_cast<int_type>(std::min(x.size(), y.size()));
          real_type * x_ = new real_type[n];
          real_type * y_ = new real_type[n];
          for (int_type i = 0; i < n; i++) {
            x_[i] = x[i];
            y_[i] = y[i];
          }
          bool ret = self->build_G1(n, x_, y_);
          delete[] x_;
          delete[] y_;
          return ret;
        })
        .def("build_G1", [](BiarcList * self, std::vector<real_type> x, std::vector<real_type> y, std::vector<real_type> theta) {
          int_type n = static_cast<int_type>(std::min({x.size(), y.size(), theta.size()}));
          real_type * x_ = new real_type[n];
          real_type * y_ = new real_type[n];
          real_type * theta_ = new real_type[n];
          for (int_type i = 0; i < n; i++) {
            x_[i] = x[i];
            y_[i] = y[i];
            theta_[i] = theta[i];
          }
          bool ret = self->build_G1(n, x_, y_);
          delete[] x_;
          delete[] y_;
          delete[] theta_;
          return ret;
        })
        .def("getAtS", &BiarcList::getAtS)
        .def("numSegments", &BiarcList::numSegments)
        .def("findAtS", &BiarcList::findAtS)
        .def("segment_length", &BiarcList::segment_length)
        .def("segment_length_ISO", &BiarcList::segment_length_ISO)
        .def("bbTriangles", [](BiarcList * self, real_type max_angle = Utils::m_pi/6, real_type max_size  = 1e100) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles(tvec, max_angle, max_size);
          return tvec;
        }, py::arg("max_angle") = Utils::m_pi/6, py::arg("max_size") = 1e100)
        .def("bbTriangles_ISO", [](BiarcList * self, real_type offs, real_type max_angle, real_type max_size) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_ISO(offs, tvec, max_angle, max_size);
          return tvec;
        })
        .def("bbTriangles_SAE", [](BiarcList * self, real_type offs, real_type max_angle = Utils::m_pi/18, real_type max_size  = 1e100) {
          std::vector<Triangle2D> tvec;
          self->bbTriangles_SAE(offs, tvec, max_angle, max_size);
          return tvec;
        }, py::arg("offs"), py::arg("max_angle") = Utils::m_pi/18, py::arg("max_size") = 1e100)
        .def("build_AABBtree_ISO", &BiarcList::build_AABBtree_ISO)
        .def("build_AABBtree_SAE", &BiarcList::build_AABBtree_SAE, 
          py::arg("offs"), py::arg("max_angle") = Utils::m_pi/6, py::arg("max_size") = 1e100)
        .def("getSTK", [](BiarcList * self) {
          int_type n = self->numSegments();
          std::vector<real_type> s, t, k;
          real_type * s_ = new real_type[n];
          real_type * t_ = new real_type[n];
          real_type * k_ = new real_type[n];
          self->getSTK(s_, t_, k_);
          for (int_type i = 0; i < n; i++) {
            s.push_back(s_[i]);
            t.push_back(t_[i]);
            k.push_back(k_[i]);
          }
          delete[] s_;
          delete[] t_;
          delete[] k_;
          return std::make_tuple(s, t, k);
        })
        .def("getXY", [](BiarcList * self) {
          int_type n = self->numSegments();
          std::vector<real_type> x, y;
          real_type * x_ = new real_type[n];
          real_type * y_ = new real_type[n];
          self->getXY(x_, y_);
          for (int_type i = 0; i < n; i++) {
            x.push_back(x_[i]);
            y.push_back(y_[i]);
          }
          delete[] x_;
          delete[] y_;
          return std::make_tuple(x, y);
        })
        .def("findST1", [](BiarcList * self, real_type x, real_type y) {
          real_type s, t;
          self->findST1(x, y, s, t);
          return std::make_tuple(s, t);
        })
        .def("findST1", [](BiarcList * self, int_type ibegin, int_type iend, real_type x, real_type y) {
          real_type s, t;
          self->findST1(ibegin, iend, x, y, s, t);
          return std::make_tuple(s, t);
        });
    }
  }
}