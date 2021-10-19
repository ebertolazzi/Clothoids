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
      py::class_<Biarc, BaseCurve>(m, "Biarc",
      R"S(
        Compute biarc fitting by Hermite data

        There are several possible constructors for this class:

         * Empty constructor
         * Copy contructor
         * Constructor from data
         * Explicit constructor from another base curve (unused)

        For the copy constructor

        :param Biarc curve: curve for copy

        For the data constructor:

        :param float x0: **x** coordinate of the origin
        :param float y0: **y** coordinate of the origin
        :param float theta0: tangent in the origin
        :param float x1: **x** coordinate of the destination
        :param float y1: **y** coordinate of the destination
        :param float theta1: tangent in the destination
      )S")

      .def(py::init())
      .def(py::init<Biarc const &>(), py::arg("curve"))
      .def(py::init<real_type, real_type, real_type, real_type, real_type, real_type>(),
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), 
        py::arg("x1"), py::arg("y1"), py::arg("theta1"))
      .def(py::init<BaseCurve const &>(), py::arg("curve"))
      
      .def("copy", &Biarc::copy, py::arg("curve"),
      R"S(
        Copy onto this object another biarc parameters

        :param Biarc curve: another biarc
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("build", &Biarc::build, 
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), 
        py::arg("x1"), py::arg("y1"), py::arg("theta1"),
      R"S(
         Construct a biarc passing from the points
         :math:`(x_0,y_0)` to the point  :math:`(x_1,y_1)`
         with initial angle :math:`\theta_0` and final angle :math:`\theta_1`
        
        :param float x0: **x** coordinate of the origin
        :param float y0: **y** coordinate of the origin
        :param float theta0: tangent in the origin
        :param float x1: **x** coordinate of the destination
        :param float y1: **y** coordinate of the destination
        :param float theta1: tangent in the destination
        :return: false if biarc cannot be computed
        :rtype: bool
      )S")

      .def("c0", &Biarc::C0, 
      R"S(
        Returns the first circle arc of the biarc

        :return: circle arc of the biarc
        :rtype: CircleArc
      )S")

      .def("c1", &Biarc::C1, 
      R"S(
        Returns the second circle arc of the biarc

        :return: circle arc of the biarc
        :rtype: CircleArc
      )S")

      .def("build_3P", &Biarc::build_3P,
        py::arg("x0"), py::arg("y0"), py::arg("x1"), py::arg("y1"), py::arg("x2"), py::arg("y2"),
      R"S(
        Construct a biarc by 3 point at "minimum energy"
        
         * Planar point set fairing and fitting by arc splines
         * Xunnian Yang and Guozhao Wang
         * Computer-Aided Design, vol 33, 2001
        
        :param float x0: **x** coordinate of initial point
        :param float y0: **y** coordinate of initial point
        :param float x1: **x** coordinate of central point
        :param float y1: **y** coordinate of central point
        :param float x2: **x** coordinate of final point
        :param float y2: **y** coordinate of final point
        :return: false if biarc cannot be computed
        :rtype: bool
      )S")
      
      .def("xMiddle", &Biarc::xMiddle, 
      R"S(
        Returns the **x** coordinate of the junction point of the biarc

        :return: **x** coordinates of the junction
        :rtype: float
      )S")
      
      .def("yMiddle", &Biarc::yMiddle,
      R"S(
        Returns the **y** coordinate of the junction point of the biarc

        :return: **y** coordinates of the junction
        :rtype: float
      )S")

      .def("thetaMiddle", &Biarc::thetaMiddle,
      R"S(
        Returns the **x** coordinate of the junction point of the biarc

        :return: **x** coordinates of the junction
        :rtype: float
      )S")
      
      .def("kappa0", &Biarc::kappa0, 
      R"S(
        Curvature of the first circle arc

        :return: curvature of the first circle arc
        :rtype: float
      )S")

      .def("kappa1", &Biarc::kappa1, 
      R"S(
        Curvature of the second circle arc

        :return: curvature of the second circle arc
        :rtype: float
      )S")

      .def("length0", &Biarc::length0, 
      R"S(
        Length of the first circle arc

        :return: length of the first circle arc
        :rtype: float
      )S")

      .def("length1", &Biarc::length1, 
      R"S(
        Length of the second circle arc

        :return: length of the second circle arc
        :rtype: float
      )S")
      
      .def("delta_theta", &Biarc::delta_theta,
      R"S(
        Change of angle of the biarc :math:`\theta_1 - \theta_0`

        :return: the angle difference
        :rtype: float
      )S");
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