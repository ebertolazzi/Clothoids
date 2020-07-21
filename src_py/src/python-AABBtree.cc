/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-AABBtree.hh"
#include <pybind11/stl.h>
#include <memory>
#include <string>
#include <sstream>

namespace G2lib {
  namespace python {
    void wrap_AABBtree(py::module & m) {
      using PyPtrBBox = std::shared_ptr<BBox>;
      using PyVecPtrBBox = std::vector<PyPtrBBox>;
      using PyPairBBox = std::pair<PyPtrBBox, PyPtrBBox>;
      using PyVecPairBBox = std::vector<PyPairBBox>;

      py::class_<BBox, PyPtrBBox>(m, "BBox")
        .def(py::init<real_type, real_type, real_type, real_type, int_type, int_type>())
        .def(py::init([](PyVecPtrBBox bboxes, int_type id, int_type ipos) {
          std::vector<BBox::PtrBBox> _bboxes;
          for (auto & el: bboxes) {
            _bboxes.push_back(static_cast<BBox::PtrBBox>(el));
          }
          return std::make_shared<BBox>(_bboxes, id, ipos);
        }))
        .def("Xmin", &BBox::Xmin)
        .def("Ymin", &BBox::Ymin)
        .def("Xmax", &BBox::Xmax)
        .def("Ymax", &BBox::Ymax)
        .def("Id", &BBox::Id)
        .def("Ipos", &BBox::Ipos)
        .def("collision", &BBox::collision)
        .def("distance", &BBox::distance)
        .def("maxDistance", &BBox::maxDistance)
        .def("join", [](PyPtrBBox self, PyVecPtrBBox bboxes) {
          std::vector<BBox::PtrBBox> _bboxes;
          for (auto & el: bboxes) {
            _bboxes.push_back(static_cast<BBox::PtrBBox>(el));
          }
          self->join(_bboxes);
        })
        .def("__str__", [](PyPtrBBox self) {
          std::ostringstream str;
          self->print(str);
          return str.str();
        });


      py::class_<AABBtree>(m, "AABBtree")
        .def(py::init<>())
        .def("clear", &AABBtree::clear)
        .def("empty", &AABBtree::empty)
        .def("bbox", [](AABBtree * self) {
          real_type x_min, y_min, x_max, y_max;
          self->bbox(x_min, y_min, x_max, y_max);
          return std::make_tuple(
            std::make_tuple(x_min, y_min), 
            std::make_tuple(x_max, y_max));
        })
        .def("build", [](AABBtree * self, PyVecPtrBBox bboxes) {
          std::vector<BBox::PtrBBox> _bboxes;
          for (auto & el: bboxes) {
            _bboxes.push_back(static_cast<BBox::PtrBBox>(el));
          }
          self->build(_bboxes);
        })
        .def("intersect", [](AABBtree * self, AABBtree const & tree, bool swap_tree) {
          AABBtree::VecPairPtrBBox intersectionList;
          PyVecPairBBox result;
          self->intersect(tree, intersectionList, swap_tree);
          for (auto & el: intersectionList) {
            BBox::PtrBBox left = std::get<0>(el);
            BBox::PtrBBox right = std::get<1>(el);
            PyPtrBBox _left = std::make_shared<BBox>(left->Xmin(), left->Ymin(), left->Xmax(), left->Ymax(), left->Id(), left->Ipos());
            PyPtrBBox _right = std::make_shared<BBox>(right->Xmin(), right->Ymin(), right->Xmax(), right->Ymax(), right->Id(), right->Ipos());
            result.push_back(std::make_pair(_left, _right));
          }
          return result;
        })
        .def("min_distance", [](AABBtree * self, real_type x, real_type y) {
          AABBtree::VecPtrBBox candidateList;
          PyVecPtrBBox result;
          self->min_distance(x, y, candidateList);
          for (auto & el: candidateList) {
            result.push_back(std::make_shared<BBox>(el->Xmin(), el->Ymin(), el->Xmax(), el->Ymax(), el->Id(), el->Ipos()));
          }
          return result;
        })
        .def("__str__", [](AABBtree * self) {
          std::ostringstream str;
          self->print(str, 0);
          return str.str();
        })
        .def("print", [](AABBtree * self, int_type level = 0) {
          std::ostringstream str;
          self->print(str, 0);
          return str.str();
        }, py::arg("level") = 0);
        // TODO: Missing collisions. Maybe I should check with Enrico
        // what I should implement as template specialization. But as for now I
        // have to wait until I implement BiarcList, ClothoidCurve, ClothoidList and Polyline
    }
  }
}