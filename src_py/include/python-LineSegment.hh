/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#pragma once

#include "pybind11/pybind11.h"
#include "python-G2libHeaders.hh"

namespace py = pybind11;

namespace G2lib {
  namespace python {
    void wrap_LineSegment(py::module & m);
    void wrap_PolyLine(py::module & m);
  }
}