/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-BaseCurve.hh"
#include "Line.hh"
#include "Circle.hh"
#include "Biarc.hh"
#include "BiarcList.hh"

namespace G2lib {
  namespace python {
    void wrap_BaseCurve(py::module & m) {
      py::enum_<CurveType>(m, "CurveType")
        .value("G2LIB_LINE", G2lib::CurveType::G2LIB_LINE)
        .value("G2LIB_POLYLINE", G2lib::CurveType::G2LIB_POLYLINE)
        .value("G2LIB_CIRCLE", G2lib::CurveType::G2LIB_CIRCLE)
        .value("G2LIB_BIARC", G2lib::CurveType::G2LIB_BIARC)
        .value("G2LIB_BIARC_LIST", G2lib::CurveType::G2LIB_BIARC_LIST)
        .value("G2LIB_CLOTHOID", G2lib::CurveType::G2LIB_CLOTHOID)
        .value("G2LIB_CLOTHOID_LIST", G2lib::CurveType::G2LIB_CLOTHOID_LIST)
        .export_values();

      py::class_<BaseCurve, PythonicBaseCurve>(m, "BaseCurve")
        .def(py::init<CurveType const &>())
        .def("type", &BaseCurve::type)
        .def("length", &BaseCurve::length)
        .def("length_ISO", &BaseCurve::length_ISO)
        .def("length_SAE", &BaseCurve::length_SAE)
        .def("thetaBegin", &BaseCurve::thetaBegin)
        .def("theta", &BaseCurve::theta)
        .def("thetaEnd", &BaseCurve::thetaEnd)
        .def("kappaBegin", &BaseCurve::kappaBegin)
        .def("kappa", &BaseCurve::kappa)
        .def("kappaEnd", &BaseCurve::kappaEnd)
        .def("xBegin", &BaseCurve::xBegin)
        .def("X", &BaseCurve::X)
        .def("xEnd", &BaseCurve::xEnd)
        .def("yBegin", &BaseCurve::yBegin)
        .def("Y", &BaseCurve::Y)
        .def("yEnd", &BaseCurve::yEnd)
        .def("xBegin_ISO", &BaseCurve::xBegin_ISO)
        .def("X_ISO", &BaseCurve::X_ISO)
        .def("xEnd_ISO", &BaseCurve::xEnd_ISO)
        .def("yBegin_ISO", &BaseCurve::yBegin_ISO)
        .def("Y_ISO", &BaseCurve::Y_ISO)
        .def("yEnd_ISO", &BaseCurve::yEnd_ISO)
        .def("xBegin_SAE", &BaseCurve::xBegin_SAE)
        .def("X_SAE", &BaseCurve::X_SAE)
        .def("xEnd_SAE", &BaseCurve::xEnd_SAE)
        .def("yBegin_SAE", &BaseCurve::yBegin_SAE)
        .def("Y_SAE", &BaseCurve::Y_SAE)
        .def("yEnd_SAE", &BaseCurve::yEnd_SAE)
        .def("tx_Begin", &BaseCurve::tx_Begin)
        .def("tx", &BaseCurve::tx)
        .def("tx_End", &BaseCurve::tx_End)
        .def("ty_Begin", &BaseCurve::ty_Begin)
        .def("ty", &BaseCurve::ty)
        .def("ty_End", &BaseCurve::ty_End)
        .def("nx_Begin_ISO", &BaseCurve::nx_Begin_ISO)
        .def("ny_Begin_ISO", &BaseCurve::ny_Begin_ISO)
        .def("nx_End_ISO", &BaseCurve::nx_End_ISO)
        .def("ny_End_ISO", &BaseCurve::ny_End_ISO)
        .def("nx_Begin_SAE", &BaseCurve::nx_Begin_SAE)
        .def("ny_Begin_SAE", &BaseCurve::ny_Begin_SAE)
        .def("nx_End_SAE", &BaseCurve::nx_End_SAE)
        .def("ny_End_SAE", &BaseCurve::ny_End_SAE)
        .def("theta_D", &BaseCurve::theta_D)
        .def("theta_DD", &BaseCurve::theta_DD)
        .def("theta_DDD", &BaseCurve::theta_DDD)
        .def("kappa_D", &BaseCurve::kappa_D)
        .def("kappa_DD", &BaseCurve::kappa_DD)
        .def("tx_D", &BaseCurve::tx_D)
        .def("tx_DD", &BaseCurve::tx_DD)
        .def("tx_DDD", &BaseCurve::tx_DDD)
        .def("ty_D", &BaseCurve::ty_D)
        .def("ty_DD", &BaseCurve::ty_DD)
        .def("ty_DDD", &BaseCurve::ty_DDD)
        .def("nx_ISO", &BaseCurve::nx_ISO)
        .def("nx_ISO_D", &BaseCurve::nx_ISO_D)
        .def("nx_ISO_DD", &BaseCurve::nx_ISO_DD)
        .def("nx_ISO_DDD", &BaseCurve::nx_ISO_DDD)
        .def("ny_ISO", &BaseCurve::ny_ISO)
        .def("ny_ISO_D", &BaseCurve::ny_ISO_D)
        .def("ny_ISO_DD", &BaseCurve::ny_ISO_DD)
        .def("ny_ISO_DDD", &BaseCurve::ny_ISO_DDD)
        .def("nx_SAE", &BaseCurve::nx_SAE)
        .def("nx_SAE_D", &BaseCurve::nx_SAE_D)
        .def("nx_SAE_DD", &BaseCurve::nx_SAE_DD)
        .def("nx_SAE_DDD", &BaseCurve::nx_SAE_DDD)
        .def("ny_SAE", &BaseCurve::ny_SAE)
        .def("ny_SAE_D", &BaseCurve::ny_SAE_D)
        .def("ny_SAE_DD", &BaseCurve::ny_SAE_DD)
        .def("ny_SAE_DDD", &BaseCurve::ny_SAE_DDD)
        .def("X_D", &BaseCurve::X_D)
        .def("X_DD", &BaseCurve::X_DD)
        .def("X_DDD", &BaseCurve::X_DDD)
        .def("Y_D", &BaseCurve::Y_D)
        .def("Y_DD", &BaseCurve::Y_DD)
        .def("Y_DDD", &BaseCurve::Y_DDD)
        .def("X_ISO_D", &BaseCurve::X_ISO_D)
        .def("X_ISO_DD", &BaseCurve::X_ISO_DD)
        .def("X_ISO_DDD", &BaseCurve::X_ISO_DDD)
        .def("Y_ISO_D", &BaseCurve::Y_ISO_D)
        .def("Y_ISO_DD", &BaseCurve::Y_ISO_DD)
        .def("Y_ISO_DDD", &BaseCurve::Y_ISO_DDD)
        .def("X_SAE_D", &BaseCurve::X_SAE_D)
        .def("X_SAE_DD", &BaseCurve::X_SAE_DD)
        .def("X_SAE_DDD", &BaseCurve::X_SAE_DDD)
        .def("Y_SAE_D", &BaseCurve::Y_SAE_D)
        .def("Y_SAE_DD", &BaseCurve::Y_SAE_DD)
        .def("Y_SAE_DDD", &BaseCurve::Y_SAE_DDD)
        .def("translate", &BaseCurve::translate)
        .def("rotate", &BaseCurve::rotate)
        .def("scale", &BaseCurve::scale)
        .def("reverse", &BaseCurve::reverse)
        .def("changeOrigin", &BaseCurve::changeOrigin)
        .def("trim", &BaseCurve::trim)
        .def("collision", &BaseCurve::collision)
        .def("collision_ISO", &BaseCurve::collision_ISO)
        .def("collision_SAE", &BaseCurve::collision_SAE)
        .def("distance", &BaseCurve::distance)
        .def("distance_ISO", &BaseCurve::distance_ISO)
        .def("distance_SAE", &BaseCurve::distance_SAE)
        .def("info", &BaseCurve::info)
        // To be modified with a lambda!
        .def("bbox", [](BaseCurve * self) {
          real_type x_min, y_min, x_max, y_max;
          self->bbox(x_min, y_min, x_max, y_max);
          return std::make_tuple(
            std::make_tuple(x_min, y_min), 
            std::make_tuple(x_max, y_max));
        })
        .def("bbox_ISO", [](BaseCurve * self, real_type offs) {
          real_type x_min, y_min, x_max, y_max;
          self->bbox_ISO(offs, x_min, y_min, x_max, y_max);
          return std::make_tuple(
            std::make_tuple(x_min, y_min), 
            std::make_tuple(x_max, y_max));
        })
        .def("bbox_SAE", [](BaseCurve * self, real_type offs) {
          real_type x_min, y_min, x_max, y_max;
          self->bbox_SAE(offs, x_min, y_min, x_max, y_max);
          return std::make_tuple(
            std::make_tuple(x_min, y_min), 
            std::make_tuple(x_max, y_max));
        })
        .def("tg", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->tg(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("tg_D", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->tg_D(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("tg_DD", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->tg_DD(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("tg_DDD", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->tg_DDD(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("nor_ISO", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->nor_ISO(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("nor_ISO_D", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->nor_ISO_D(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("nor_ISO_DD", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->nor_ISO_DD(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("nor_ISO_DDD", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->nor_ISO_DDD(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("nor_SAE", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->nor_SAE(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("nor_SAE_D", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->nor_SAE_D(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("nor_SAE_DD", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->nor_SAE_DD(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("nor_SAE_DDD", [](BaseCurve * self, real_type s) {
          real_type r_x, r_y;
          self->nor_SAE_DDD(s, r_x, r_y);
          return std::make_tuple(r_x, r_y);
        })
        .def("evaluate", [](BaseCurve * self, real_type s) {
          real_type th, k, x, y;
          self->evaluate(s, th, k, x, y);
          return std::make_tuple(th, k, x, y);
        })
        .def("eval", [](BaseCurve * self, real_type s) {
          real_type x, y;
          self->eval(s, x, y);
          return std::make_tuple(x, y);
        })
        .def("evaluate_ISO", [](BaseCurve * self, real_type s, real_type offs) {
          real_type th, k, x, y;
          self->evaluate_ISO(s, offs, th, k, x, y);
          return std::make_tuple(th, k, x, y);
        })
        .def("eval_ISO", [](BaseCurve * self, real_type s, real_type offs) {
          real_type x, y;
          self->eval_ISO(s, offs, x, y);
          return std::make_tuple(x, y);
        })
        .def("evaluate_SAE", [](BaseCurve * self, real_type s, real_type offs) {
          real_type th, k, x, y;
          self->evaluate_SAE(s, offs, th, k, x, y);
          return std::make_tuple(th, k, x, y);
        })
        .def("eval_SAE", [](BaseCurve * self, real_type s, real_type offs) {
          real_type x, y;
          self->eval_SAE(s, offs, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_D", [](BaseCurve * self, real_type s) {
          real_type x, y;
          self->eval_D(s, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_DD", [](BaseCurve * self, real_type s) {
          real_type x, y;
          self->eval_DD(s, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_DDD", [](BaseCurve * self, real_type s) {
          real_type x, y;
          self->eval_DDD(s, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_ISO_D", [](BaseCurve * self, real_type s, real_type offs) {
          real_type x, y;
          self->eval_ISO_D(s, offs, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_ISO_DD", [](BaseCurve * self, real_type s, real_type offs) {
          real_type x, y;
          self->eval_ISO_DD(s, offs, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_ISO_DDD", [](BaseCurve * self, real_type s, real_type offs) {
          real_type x, y;
          self->eval_ISO_DDD(s, offs, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_SAE_D", [](BaseCurve * self, real_type s, real_type offs) {
          real_type x, y;
          self->eval_SAE_D(s, offs, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_SAE_DD", [](BaseCurve * self, real_type s, real_type offs) {
          real_type x, y;
          self->eval_SAE_DD(s, offs, x, y);
          return std::make_tuple(x, y);
        })
        .def("eval_SAE_DDD", [](BaseCurve * self, real_type s, real_type offs) {
          real_type x, y;
          self->eval_SAE_DDD(s, offs, x, y);
          return std::make_tuple(x, y);
        })
        .def("findST_ISO", [](BaseCurve * self, real_type x, real_type y) {
          real_type s, t;
          self->findST_ISO(x, y, s, t);
          return std::make_tuple(s, t);
        })
        .def("findST_SAE", [](BaseCurve * self, real_type x, real_type y) {
          real_type s, t;
          self->findST_SAE(x, y, s, t);
          return std::make_tuple(s, t);
        })
        .def("closestPoint_ISO", [](BaseCurve * self, real_type qx, real_type qy) {
          real_type x, y, s, t, dst;
          self->closestPoint_ISO(qx, qy, x, y, s, t, dst);
          return std::make_tuple(x, y, s, t, dst);
        })  
        .def("closestPoint_ISO", [](BaseCurve * self, real_type qx, real_type qy, real_type offs) {
          real_type x, y, s, t, dst;
          self->closestPoint_ISO(qx, qy, offs, x, y, s, t, dst);
          return std::make_tuple(x, y, s, t, dst);
        })   
        .def("closestPoint_SAE", [](BaseCurve * self, real_type qx, real_type qy) {
          real_type x, y, s, t, dst;
          self->closestPoint_SAE(qx, qy, x, y, s, t, dst);
          return std::make_tuple(x, y, s, t, dst);
        })  
        .def("closestPoint_SAE", [](BaseCurve * self, real_type qx, real_type qy, real_type offs) {
          real_type x, y, s, t, dst;
          self->closestPoint_SAE(qx, qy, offs, x, y, s, t, dst);
          return std::make_tuple(x, y, s, t, dst);
        })
        .def("__str__", [](BaseCurve * self) {
          std::ostringstream str;
          self->info(str);
          return str.str();
        });
      
      // py::class_<G2lib::Biarc, G2lib::BaseCurve>(m, "Biarc")
      //   .def(py::init<>())
      //   .def(py::init<G2lib::BaseCurve const &>())
      //   .def(py::init<G2lib::Biarc const &>())
      //   .def(py::init<real_type, real_type, real_type, real_type, real_type, real_type>());

      // py::class_<G2lib::BiarcList, G2lib::BaseCurve>(m, "BiarcList")
      //   .def(py::init<>())
      //   .def(py::init<G2lib::BaseCurve const &>())
      //   .def(py::init<G2lib::BiarcList const &>())
      //   .def(py::init<G2lib::LineSegment const &>())
      //   .def(py::init<G2lib::CircleArc const &>());     
    }
  }
}