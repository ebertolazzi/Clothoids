/**
 * PYTHON Wrapper for Clothoids
 * 
 * License MIT - See LICENSE file
 * 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
 *      Enrico Bertolazzi, Marco Frego
 */

#include "python-Clothoid.hh"
#include "python-ClothoidSpline-Interpolation.hh"
#include "pybind11/stl.h"
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
      py::class_<ClothoidCurve, BaseCurve>(m, "ClothoidCurve",
      R"S(
        Class that menages a clothoid curve (a curve with linearly variable 
        curvature). There are several possible constructors for this class:

         * constructor from a base curve
         * constructor from a line segment
         * constructor from a circle arc
         * constructor from data
         * G1 constructor

        For the constructor from data it is possible to use the following
        parameters:

        :param float x0: starting position **x** coordinate
        :param float y0: starting position **y** coordinate
        :param float theta0: initial angle
        :param float k: initial curvature
        :param float dk: curvature derivative
        :param float L: length of the curve

        Instead, if one wants to create the clothoid by solving the G1 problem:

        :param Tuple[float, float] p0: starting position **(x, y)** coordinates
        :param float theta0: initial angle
        :param Tuple[float, float] p1: ending position **(x, y)** coordinates
        :param float theta1: final angle
      )S")
        
      .def(py::init())
      .def(py::init<BaseCurve const &>())
      .def(py::init<ClothoidCurve const &>())
      .def(py::init<LineSegment const &>())
      .def(py::init<CircleArc const &>())
      .def(py::init<real_type, real_type, real_type, real_type, real_type, real_type>(),
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("k"), py::arg("dk"), py::arg("L"))
      .def(py::init([](std::tuple<real_type, real_type> P0, real_type theta0, 
                       std::tuple<real_type, real_type> P1, real_type theta1){
        real_type _P0[2] = {std::get<0>(P0), std::get<1>(P0)};
        real_type _P1[2] = {std::get<0>(P1), std::get<1>(P1)};
        return ClothoidCurve(_P0, theta0, _P1, theta1);
      }), py::arg("p0"), py::arg("theta0"), py::arg("p1"), py::arg("theta1"))

      .def("copy", [](const ClothoidCurve & self) {
        ClothoidCurve other;
        other.copy(self);
        return other;
      },
      R"S(
        Create a copy of the current clothoid curve

        :return: a new copy of the current clothoid
        :rtype: ClothoidCurve
      )S")

      .def("build_G1", &ClothoidCurve::build_G1, 
        py::arg("x0"), py::arg("y0"), py::arg("theta0"),
        py::arg("x1"), py::arg("y1"), py::arg("theta1"), py::arg("tol") = 1e-12,
      R"S(
        Build a clothoid by solving the hermite G1 problem
        
        :param float x0:     initial x position :math:`x_0`
        :param float y0:     initial y position :math:`y_0`
        :param float theta0: initial angle      :math:`\theta_0`
        :param float x1:     final x position   :math:`x_1`
        :param float y1:     final y position   :math:`y_1`
        :param float theta1: initial angle      :math:`\theta_1`
        :param float tol: tolerance (default 1e-12)
        :return: number of iteration performed
        :rtype: int
      )S")
      
      .def("build", py::overload_cast<real_type, real_type, real_type, real_type, real_type, real_type>(&ClothoidCurve::build),
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("k"), py::arg("dk"), py::arg("L"),
      R"S(
        Builds a clothoid curve starting from raw data

        :param float x0: starting position **x** coordinate
        :param float y0: starting position **y** coordinate
        :param float theta0: initial angle
        :param float k: initial curvature
        :param float dk: curvature derivative
        :param float L: length of the curve
        :return: nothing, works in place 
        :rtype: NoneType
      )S")

      .def("build_G1_D", [](ClothoidCurve * self,
                            real_type x0, real_type y0, real_type theta0,
                            real_type x1, real_type y1, real_type theta1,
                            std::tuple<real_type, real_type> L_D, std::tuple<real_type, real_type> k_D, 
                            std::tuple<real_type, real_type> dk_D, real_type tol = 1e-12) {
        real_type _L_D[2] = {std::get<0>(L_D), std::get<1>(L_D)};
        real_type _k_D[2] = {std::get<0>(k_D), std::get<1>(k_D)};
        real_type _dk_D[2] = {std::get<0>(dk_D), std::get<1>(dk_D)};
        return self->build_G1_D(x0, y0, theta0, x1, y1, theta1, _L_D, _k_D, _dk_D, tol);
      }, 
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), 
        py::arg("x1"), py::arg("y1"), py::arg("theta1"),
        py::arg("L_D"), py::arg("k_D"), py::arg("dk_D"), py::arg("tol") = 1e-12,
      R"S(
        Build a clothoid by solving the hermite G1 problem.
        
        :param float x0:     initial x position :math:`x_0`
        :param float y0:     initial y position :math:`y_0`
        :param float theta0: initial angle      :math:`\theta_0`
        :param float x1:     final x position   :math:`x_1`
        :param float y1:     final y position   :math:`y_1`
        :param float theta1: final angle        :math:`\theta_1`
        :param Tuple[float, float] L_D:    derivative of the length :math:`L(\theta_0,\theta_1)`
        :param Tuple[float, float] k_D:    derivative of the curvature :math:`\kappa(\theta_0,\theta_1)`
        :param Tuple[float, float] dk_D:   derivative of the curvature variation :math:`\kappa'(\theta_0,\theta_1)`
        :param float tol: tolerance (default 1e-12)
        :return: number of iteration performed
        :rtype: int
      )S")

      .def("build_forward", &ClothoidCurve::build_forward, 
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("kappa0"),
        py::arg("x1"), py::arg("y1"), py::arg("tol") = 1e-12,
      R"S(
        Build a clothoid by solving the forward problem
        
        :param float x0:     initial x position :math:`x_0`
        :param float y0:     initial y position :math:`y_0`
        :param float theta0: initial angle      :math:`\theta_0`
        :param float x1:     final x position   :math:`x_1`
        :param float y1:     final y position   :math:`y_1`
        :param float tol: tolerance (default 1e-12)
        :return: true if solution exists
        :rtype: bool
      )S")

      .def("build", py::overload_cast<LineSegment const &>(&ClothoidCurve::build),
        py::arg("ls"),
      R"S(
        Builds a clothoid curve starting from a line segment

        :param LineSegment ls: the line segment
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("build", py::overload_cast<CircleArc const &>(&ClothoidCurve::build),
        py::arg("ca"),
      R"S(
        Builds a clothoid curve starting from a circle arc

        :param CircleArc ls: the circle arc
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("Pinfinity", [](const ClothoidCurve & self, bool plus) {
        real_type x, y;
        self.Pinfinity(x, y, plus);
        return std::make_tuple(x, y);
      }, py::arg("plus") = true,
      R"S(
        Return the point at infinity of the clothoids :math:`P(s)`.
    
        :param float x:    x-coordinate of the infinity point
        :param float y:    y-coordinate of the infinity point
        :param float plus: it true return :math:`\lim_{s\to+\infty} P(s) `
                      otherwise return :math:`\lim_{s\to-\infty} P(s)`
        :return: the point coordinates as a tuple
        :rtype: Tuple[float, float]
      )S")
      
      .def("dkappa", &ClothoidCurve::dkappa,
      R"S(
        Derivative of the curvature of the clothoid

        :return: derivative of the curvature
        :rtype: float
      )S")
        
      .def("thetaTotalVariation", &ClothoidCurve::thetaTotalVariation,
      R"S(
        Clothoid curve total angle variation

        :return: total angle variation
        :rtype: float
      )S")

      .def("thetaMinMax", [](const ClothoidCurve & self) {
        real_type th_min, th_max;
        self.thetaMinMax(th_min, th_max);
        return std::make_tuple(th_min, th_max);
      },
      R"S(
        Returns the theta minimum and maximum for the curve

        :return: theta minimum and maximum
        :rtype: Tuple[float, float]
      )S")

      .def("deltaTheta", &ClothoidCurve::deltaTheta,
      R"S(
        Returns the difference between theta minimum and theta maximum

        :return: :math:`\theta_{\max} - \thate_{\min}`
        :rtype: float
      )S")

      .def("curvatureMinMax", [](const ClothoidCurve & self) {
        real_type cv_min, cv_max;
        self.curvatureMinMax(cv_min, cv_max);
        return std::make_tuple(cv_min, cv_max);
      }, 
      R"S(
        Curvature minimum and maximum value for the curve

        :return: minimum and maximum in curvature
        :rtype: Tuple[float, float]
      )S")

      .def("curvatureTotalVariation", &ClothoidCurve::curvatureTotalVariation,
      R"S(
        Clothoid total curvature variation

        :return: clothoid total curvature variation
        :rtype: float
      )S")
      
      .def("integralCurvature2", &ClothoidCurve::integralCurvature2,
      R"S(
        Given the clothoid curve :math:`P(s)` compute.
     
        .. math: \int_0^L |P''(s)|^2 \mathrm{d}s

        :return: integral curvature
        :rtype: float
      )S")
      
      .def("integralJerk2", &ClothoidCurve::integralJerk2,
      R"S(
        Given the clothoid curve :math:`P(s)` compute.
     
        .. math: \int_0^L |P'''(s)|^2 \mathrm{d}s

        :return: integral jerk
        :rtype: float
      )S")
      
      .def("integralSnap2", &ClothoidCurve::integralSnap2,
      R"S(
        Given the clothoid curve :math:`P(s)` compute.
     
        .. math: \int_0^L |P''''(s)|^2 \mathrm{d}s

        :return: integral jerk
        :rtype: float
      )S")
      
      .def("optimized_sample_ISO", [](const ClothoidCurve & self, real_type offs, int_type npts, real_type max_angle) {
        std::vector<real_type> ret;
        self.optimized_sample_ISO(offs, npts, max_angle, ret);
        return ret;
      }, py::arg("offs"), py::arg("npts"), py::arg("max_angle"),
      R"S(
        Return a vector of optimized sample parameters for plotting. It uses
        the ISO standard reference frame
        
        :param offs: offset of the sampled curve
        :param npts: suggested minimum number of sampled points
        :param max_angle: maximum angle variation between two sampled points
        :return: vector of computed parameters
        :rtype: List[float]
      )S")
      
      .def("optimized_sample_SAE", [](const ClothoidCurve & self, real_type offs, int_type npts, real_type max_angle) {
        std::vector<real_type> ret;
        self.optimized_sample_SAE(offs, npts, max_angle, ret);
        return ret;
      }, py::arg("offs"), py::arg("npts"), py::arg("max_angle"),
      R"S(
        Return a vector of optimized sample parameters for plotting. It uses
        the SAE standard reference frame
        
        :param offs: offset of the sampled curve
        :param npts: suggested minimum number of sampled points
        :param max_angle: maximum angle variation between two sampled points
        :return: vector of computed parameters
        :rtype: List[float]
      )S")
      
      .def("closestPointBySample", [](ClothoidCurve * self, real_type ds, int_type qx, real_type qy) {
        real_type x, y, s;
        real_type v = self->closestPointBySample(ds, qx, qy, x, y, s);
        return std::make_tuple(v, x, y, s);
      }, py::arg("ds"), py::arg("qx"), py::arg("qy"),
      R"S(
        Compute the point on clothoid at minimal distance from a given point
        using the optimized algorithm described in the publication:
        
         * **E.Bertolazzi, M.Frego**, Point-Clothoid distance and projection computation
           SIAM J. Scientific Computing, Vol. 41, No. 5, pp. A3326-A3353

        Returns a tuple containing:

         * the distance of the point from the clothoid
         * x-coordinate of the point on clothoid at minimal distance
         * y-coordinate of the point on clothoid at minimal distance
         * curvilinear coordinate of the point (X,Y) on the clothoid
        
        :param float ds: sampling step
        :param float qx: x-coordinate of the given point
        :param float qy: y-coordinate of the given point
        :return the: the tuple as described
        :rtype: Tuple[float, float, float, float]
      )S")
      
      .def("distanceBySample", [](const ClothoidCurve & self, real_type ds, int_type qx, real_type qy) {
        real_type s;
        real_type v = self.distanceBySample(ds, qx, qy, s);
        return std::make_tuple(v, s);
      }, py::arg("ds"), py::arg("qx"), py::arg("qy"),
      R"S(
        Approximate the point on clothoid at minimal distance from a given point
        using simple sampling.
         
        :param  ds: sampling step
        :param  qx: x-coordinate of the given point
        :param  qy: y-coordinate of the given point
        :return: the distance of the point from the clothoid and the 
                 curvilinear coordinate of the point (X,Y) on the clothoid
        :rtype: Tuple[float, float]
      )S")
      
      .def("changeCurvilinearOrigin", &ClothoidCurve::changeCurvilinearOrigin,
      R"S(
        Change the origin of the clothoid at :math:`s_0`
        and the length to  :math:`L`.
         
        :param float s0:   :math:`s_0`
        :param float newL: :math:`L`
        :return: nothing, works in place
        :rtype: NoneType
      )S")
      
      .def("build_AABBtree_ISO", &ClothoidCurve::build_AABBtree_ISO,
        py::arg("offs"), py::arg("max_angle"), py::arg("max_size"),
      R"S(
        Builds the AABB tree of the current curve. Uses ISO reference 
        frame for offsets

        :param float offs: offset from the curve
        :param float max_angle: maximum angle
        :param float max_size: maximum size
        :return: nothing works in place
        :rtype: NoneType
      )S")
      
      .def("approximate_collision_ISO", &ClothoidCurve::approximate_collision_ISO,
        py::arg("offs"), py::arg("c"), py::arg("offs_C"), py::arg("max_angle"), py::arg("max_size"),
      R"S(
        Collision detection. Uses ISO standard reference frame for offsets.
         
        :param float offs:      curve offset
        :param ClothoidCurve c: curve to compare for collision detection
        :param float offs_C:    curve offset
        :param float max_angle: maximum angle variation
        :param float max_size:  if the segment is larger then this parameter is split
         
      )S");
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

      m.def("buildP1", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys, real_type theta_0, real_type theta_1) -> G2lib::ClothoidList {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator interpolator(xs, ys);
        interpolator.buildP1(theta_0, theta_1, result);
        return result;
      }, py::arg("xs"), py::arg("ys"), py::arg("theta0"), py::arg("theta1"),
      R"S(
        Builds a clothoid list starting from a list of points. Build a 
        clothoids between each point pair.

        Uses target P1. Requires Eigen library during compilation

        :param List[float] xs: **x** coordinates of points
        :param List[float] ys: **y** coordinates of points
        :param float theta0: intial angle
        :param float theta1: final angle 
        :return: the clothoid list
        :rtype: ClothodList 
      )S")
      .def("buildP2", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> G2lib::ClothoidList {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator interpolator(xs, ys);
        interpolator.buildP2(result);
        return result;
      }, py::arg("xs"), py::arg("ys"), 
      R"S(
        Builds a clothoid list starting from a list of points. Build a 
        clothoids between each point pair.

        Uses target P2. Requires Eigen library during compilation

        :param List[float] xs: **x** coordinates of points
        :param List[float] ys: **y** coordinates of points
        :return: the clothoid list
        :rtype: ClothodList 
      )S")
      .def("buildP4", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> G2lib::ClothoidList {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator interpolator(xs, ys);
        interpolator.buildP4(result);
        return result;
      }, py::arg("xs"), py::arg("ys"), 
      R"S(
        Builds a clothoid list starting from a list of points. Build a 
        clothoids between each point pair.

        Uses target P4. Requires IPOPT library during compilation

        :param List[float] xs: **x** coordinates of points
        :param List[float] ys: **y** coordinates of points
        :return: the clothoid list
        :rtype: ClothodList 
      )S")
      .def("buildP5", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> G2lib::ClothoidList {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator interpolator(xs, ys);
        interpolator.buildP5(result);
        return result;
      }, py::arg("xs"), py::arg("ys"), 
      R"S(
        Builds a clothoid list starting from a list of points. Build a 
        clothoids between each point pair.

        Uses target P5. Requires IPOPT library during compilation

        :param List[float] xs: **x** coordinates of points
        :param List[float] ys: **y** coordinates of points
        :return: the clothoid list
        :rtype: ClothodList 
      )S")
      .def("buildP6", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> G2lib::ClothoidList {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator interpolator(xs, ys);
        interpolator.buildP6(result);
        return result;
      }, py::arg("xs"), py::arg("ys"), 
      R"S(
        Builds a clothoid list starting from a list of points. Build a 
        clothoids between each point pair.

        Uses target P6. Requires IPOPT library during compilation

        :param List[float] xs: **x** coordinates of points
        :param List[float] ys: **y** coordinates of points
        :return: the clothoid list
        :rtype: ClothodList 
      )S")
      .def("buildP7", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> G2lib::ClothoidList {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator interpolator(xs, ys);
        interpolator.buildP7(result);
        return result;
      }, py::arg("xs"), py::arg("ys"), 
      R"S(
        Builds a clothoid list starting from a list of points. Build a 
        clothoids between each point pair.

        Uses target P7. Requires IPOPT library during compilation

        :param List[float] xs: **x** coordinates of points
        :param List[float] ys: **y** coordinates of points
        :return: the clothoid list
        :rtype: ClothodList 
      )S")
      .def("buildP8", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> G2lib::ClothoidList {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator interpolator(xs, ys);
        interpolator.buildP8(result);
        return result;
      }, py::arg("xs"), py::arg("ys"), 
      R"S(
        Builds a clothoid list starting from a list of points. Build a 
        clothoids between each point pair.

        Uses target P8. Requires IPOPT library during compilation

        :param List[float] xs: **x** coordinates of points
        :param List[float] ys: **y** coordinates of points
        :return: the clothoid list
        :rtype: ClothodList 
      )S")
      .def("buildP9", [](const std::vector<real_type> & xs, const std::vector<real_type> & ys) -> G2lib::ClothoidList {
        G2lib::ClothoidList result;
        G2lib::Interpolation::Interpolator interpolator(xs, ys);
        interpolator.buildP9(result);
        return result;
      }, py::arg("xs"), py::arg("ys"), 
      R"S(
        Builds a clothoid list starting from a list of points. Build a 
        clothoids between each point pair.

        Uses target P9. Requires IPOPT library during compilation

        :param List[float] xs: **x** coordinates of points
        :param List[float] ys: **y** coordinates of points
        :return: the clothoid list
        :rtype: ClothodList 
      )S");
    }

    void wrap_ClothoidList(py::module & m) {
      py::class_<ClothoidList, BaseCurve>(m, "ClothoidList",
      R"S(
        Manages a piecewise clothoid composed by n clothoids, not
        necessarily G2 or G1 connected

        There are several constructors for this class:

         * empty constructor
         * copy constructor from another clothoid list
         * constructor from another base curve, such as line segment,
           circle arc, biarc, ciarc list, clothoid curve or polyline

        :param BaseCurve c: another curve for clothoid list initialization
                            (optional parameter)
      )S")

      .def(py::init<>())
      .def(py::init<ClothoidList const &>())
      .def(py::init<LineSegment const &>())
      .def(py::init<CircleArc const &>())
      .def(py::init<Biarc const &>())
      .def(py::init<BiarcList const &>())
      .def(py::init<ClothoidCurve const &>())
      .def(py::init<PolyLine const &>())
      .def(py::init<BaseCurve const &>())
      
      .def("get", &ClothoidList::get, py::arg("idx"),
      R"S(
        Returns the `idx`-th element of the list

        :param int idx: index
        :return: the curve at that index
        :rtype: BaseCurve
      )S")

      .def("__getitem__", &ClothoidList::get, py::arg("idx"),
      R"S(
        Returns the `idx`-th element of the list

        :param int idx: index
        :return: the curve at that index
        :rtype: BaseCurve
      )S")

      .def("copy", [](const ClothoidList & self) {
        ClothoidList other;
        other.copy(self);
        return other;
      }, 
      R"S(
        Returns a copy of the current biarc list

        :return: a copy of the current biarc list
        :rtype: ClothoidList
      )S")

      .def("init", &ClothoidList::init,
      R"S(
        Empties the current clothoid list

        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("reserve", &ClothoidList::reserve, py::arg("n"),
      R"S(
        Reserve memory for **n** curves in the list

        :param int n: number of element to reserve
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<LineSegment const &>(&ClothoidList::push_back), py::arg("c"),
      R"S(
        Append a line segment to the current list. The curve is transformed
        to a degenerate clothoid curve

        :param LineSegment c: the line segment
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<CircleArc const &>(&ClothoidList::push_back), py::arg("c"),
      R"S(
        Append a circle arc to the current list. The curve is transformed
        to a degenerate clothoid curve

        :param CircleArc c: the circle arc
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<Biarc const &>(&ClothoidList::push_back), py::arg("c"),
      R"S(
        Append a biarc to the current list. The curve is transformed to 
        a degenerate clothoid

        :param Biarc c: the biarc
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<BiarcList const &>(&ClothoidList::push_back), py::arg("c"),
      R"S(
        Append a biarc to the current list. The curve is transformed to 
        a degenerate clothoid list

        :param BiarcList c: the biarc list
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<PolyLine const &>(&ClothoidList::push_back), py::arg("c"),
      R"S(
        Append a polyline to the current list. The curve is transformed
        to a degenerate clothoid list

        :param PolyLine c: the polyline
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<ClothoidCurve const &>(&ClothoidList::push_back), py::arg("c"),
      R"S(
        Append a clothoid curve to the current list.

        :param ClothoidCurve c: the clothoid
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<ClothoidList const &>(&ClothoidList::push_back), py::arg("c"),
      R"S(
        Concatenate two clothoid list onto the current list.

        :param ClothoidList c: the clothoid list
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<real_type, real_type, real_type>(&ClothoidList::push_back),
        py::arg("kappa0"), py::arg("dkappa"), py::arg("l"),
      R"S(
        Adds a clothoid to the tail of the clothoid list

        :param float kappa0: initial curvature
        :param float dkappa: derivative of the curvature
        :param float l: length of the appended curve
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back", py::overload_cast<real_type, real_type, real_type, real_type, real_type, real_type>(&ClothoidList::push_back),
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("kappa0"), py::arg("dkappa"), py::arg("l"),
      R"S(
        Adds a clothoid to the tail of the clothoid list
        The builded clothoid is translated to the tail of the clothioid list.

        :param float x0: **x** coordinate of the origin
        :param float y0: **y** coordinate of the origin
        :param float theta0: tangent in the origin
        :param float kappa0: initial curvature
        :param float dkappa: derivative of the curvature
        :param float l: length of the appended curve
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back_G1", py::overload_cast<real_type, real_type, real_type>(&ClothoidList::push_back_G1),
        py::arg("x1"), py::arg("y1"), py::arg("theta1"),
      R"S(
        Add a clothoid to the tail of the clothoid list solving the G1 problem.
        The initial point and angle are taken from the tail of the clothoid list.
      
        :param flaot x1:     final x
        :param flaot y1:     final y
        :param flaot theta1: final angle
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("push_back_G1", py::overload_cast<real_type, real_type, real_type, real_type, real_type, real_type>(&ClothoidList::push_back_G1),
        py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("x1"), py::arg("y1"), py::arg("theta1"),
      R"S(
        Add a clothoid to the tail of the clothoid list solving the G1 problem.
        The initial point and angle are taken from the tail of the clothoid list.
        The builded clothoid is translated to the tail of the clothioid list.
      
        :param flaot x0:     initial x
        :param flaot y0:     initial y
        :param flaot theta0: initial angle
        :param flaot x1:     final x
        :param flaot y1:     final y
        :param flaot theta1: final angle
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("is_closed", &ClothoidList::is_closed,
      R"S(
        True if curve is closed

        :return: true if the curve is closed
        :rtype: bool
      )S")

      .def("make_closed", &ClothoidList::make_closed,
      R"S(
        Set closure flag to true

        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("make_open", &ClothoidList::make_open,
      R"S(
        Set closure flag to false

        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("closure_gap_x", &ClothoidList::closure_gap_x,
      R"S(
        Difference between initial and final points along **x**

        :return: the distance between the initial and final point along **x**
        :rtype: float
      )S")

      .def("closure_gap_y", &ClothoidList::closure_gap_y,
      R"S(
        Difference between initial and final points along **y**

        :return: the distance between the initial and final point along **y**
        :rtype: float
      )S")

      .def("closure_gap_tx", &ClothoidList::closure_gap_tx,
      R"S(
        Difference between initial and final tangents along **x**

        :return: the distance between the initial and final tangents along **x**
        :rtype: float
      )S")

      .def("closure_gap_ty", &ClothoidList::closure_gap_ty,
      R"S(
        Difference between initial and final tangents along **y**

        :return: the distance between the initial and final tangents along **y**
        :rtype: float
      )S")

      .def("closure_check", &ClothoidList::closure_check,
        py::arg("tol_xy") = 1e-6, py::arg("tol_tg") = 1e-6,
      R"S(
        Check if curve is closed (and tanget continue)

        :param float tol_xy: position tolerance
        :param float tol_tg: tangent tolerance
        :return: true if the curve is closed
        :rtype: bool
      )S")

      .def("build_G1", [](ClothoidList & self, const std::vector<real_type> & x, const std::vector<real_type> & y) {
        const size_t n = std::min(x.size(), y.size());
        return self.build_G1(n, x.data(), y.data());
      }, py::arg("x"), py::arg("y"),
      R"S(
        Build clothoid list passing to a list of points solving a series of 
        G1 fitting problems. The angle at points are estimated using the 
        routine `xy_to_guess_angle`
        
        :param List[float] x: x-coordinates list
        :param List[float] y: y-coordinates list
        :return: false if routine fails
        :rtype: bool
      )S")

      .def("build_G1", [](ClothoidList & self, const std::vector<real_type> & x, const std::vector<real_type> & y, const std::vector<real_type> & theta) {
        const size_t n = std::min({x.size(), y.size(), theta.size()});
        return self.build_G1(n, x.data(), y.data(), theta.data());
      }, py::arg("x"), py::arg("y"), py::arg("theta"),
      R"S(
        Build clothoid list passing to a list of points solving a series of 
        G1 fitting problems. Requires angles.
        
        :param List[float] x: x-coordinates list
        :param List[float] y: y-coordinates list
        :param List[float] theta: angle values
        :return: false if routine fails
        :rtype: bool
      )S")

      .def("build", [](ClothoidList & self, real_type x0, real_type y0, real_type theta0, const std::vector<real_type> & s, const std::vector<real_type> & kappa) {
        const size_t n = std::min(s.size(), kappa.size());
        return self.build(x0, y0, theta0, n, s.data(), kappa.data());
      }, py::arg("x0"), py::arg("y0"), py::arg("theta0"), py::arg("s"), py::arg("kappa"),
      R"S(
        Build clothoid list with G2 continuity.
        The vector `s` contains the breakpoints of the curve.
        Between two breakpoint the curvature change linearly (is a clothoid)
        
        :param float x0:          initial x
        :param float y0:          initial y
        :param float theta0:      initial angle
        :param List[float] s:     break point of the piecewise curve
        :param List[float] kappa: curvature at the break point
        :return: true if curve is closed
        :rtype: bool
      )S")

      .def("build_raw", [](ClothoidList & self, const std::vector<real_type> & x, const std::vector<real_type> & y,
        const std::vector<real_type> & s, const std::vector<real_type> & theta, const std::vector<real_type> & kappa) {
          return self.build_raw(x, y, s, theta, kappa);
        }, py::arg("x"), py::arg("y"), py::arg("s"), py::arg("theta"), py::arg("kappa"),
      R"S(
        Builds a clothoid list using raw data
        
        :param List[float] x:     x-coordinates
        :param List[float] y:     y-coordinates
        :param List[float] s:     break point of the piecewise curve
        :param List[float] theta: angle at break point
        :param List[float] kappa: curvature at the break point
        :return: false if fails
        :rtype: bool
      )S")

      .def("getAtS", &ClothoidList::getAtS, py::arg("s"),
      R"S(
        Get the `idx`-th clothoid of the list where `idx` is the clothoid at parameter `s`

        :param float s: curvilinear abscissa
        :return: the index of the curve at that abscissa
        :rtype: int
      )S")

      .def("numSegments", &ClothoidList::numSegments,
      R"S(
        Returns the number of segments in the list

        :return: the number of segments
        :rtype: int
      )S")

      .def("__len__", &ClothoidList::numSegments,
      R"S(
        Returns the number of segments in the list

        :return: the number of segments
        :rtype: int
      )S")

      .def("wrap_in_range", [](const ClothoidList & self, real_type s) {
        self.wrap_in_range(s);
        return s;
      }, py::arg("s"),
      R"S(
        The list of clothoid has total length :math:`L`
        the parameter :math:`s` is recomputed as :math:`s+kL` in such a way
        :math:`s+kL\in[0,L)` with :math:`k\in\mathbb{Z}`

        :param float s: the wrapping input s
        :return: the wrapped s
        :rtype: float 
      )S")

      .def("findAtS", &ClothoidList::findAtS, py::arg("s"),
      R"S(
        Find the clothoid segment whose definition range contains `s`

        :param float s: curvilinear abscissa
        :return: the segment index
        :rtype: int
      )S")

      .def("segment_length", &ClothoidList::segment_length, py::arg("nseg"),
      R"S(
        Returns the length of the segment at the requested index

        :param int nseg: index of the segment to measure
        :return: the length of the segment
        :rtype: float
      )S")

      .def("segment_length_ISO", &ClothoidList::segment_length_ISO, py::arg("nseg"), py::arg("offs"),
      R"S(
        Returns the length of the segment at the requested index, with
        an offset from the segment in ISO reference frame

        :param int nseg: index of the segment to measure
        :param float offs: offset from the segmente
        :return: the length of the segment
        :rtype: float
      )S")

      .def("segment_length_SAE", &ClothoidList::segment_length_SAE, py::arg("nseg"), py::arg("offs"),
      R"S(
        Returns the length of the segment at the requested index, with
        an offset from the segment in SAE reference frame

        :param int nseg: index of the segment to measure
        :param float offs: offset from the segmente
        :return: the length of the segment
        :rtype: float
      )S")

      .def("closestSegment", &ClothoidList::closestSegment, py::arg("qx"), py::arg("qy"),
      R"S(
        Given a point, returns the index of the closest segment to that point

        :param float qx: **x** coordinates of the point
        :param float qy: **y** coordinates of the point
        :return: the index of the closest segment
        :rtype: int
      )S")

      .def("closestPointInRange", [](const ClothoidList & self, real_type qx, real_type qy, int_type icurve_begin, int_type icurve_end) {
        real_type x, y, s, t, dst;
        int_type icurve;
        int_type ret = self.closestPointInRange_ISO(qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve);
        return std::make_tuple(ret, x, y, s, t, dst, icurve);
      }, 
        py::arg("qx"), py::arg("qy"), py::arg("icurve_begin"), py::arg("icurve_end"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the set of curve
        inside the range of segment ``icurve_begin..icurve_end``.

        There are 7 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve
         7. **index** of the closest segment 

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :param int icurve_begin: first segment of the range
        :param int icurve_end: last segment of the range
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, float, int]
      )S")

      .def("closestPointInRange", [](const ClothoidList & self, const std::vector<real_type> & qx, 
        const std::vector<real_type> & qy, int_type icurve_begin, int_type icurve_end) {
        const size_t n = std::min(qx.size(), qy.size());
        std::vector<std::tuple<int_type, real_type, real_type, real_type, real_type, real_type, int_type>> data(n);
        for (size_t i = 0; i < n; i++) {
          real_type x, y, s, t, dst;
          int_type icurve;
          int_type ret = self.closestPointInRange_ISO(qx[i], qy[i], icurve_begin, icurve_end, x, y, s, t, dst, icurve);
          data[i] = std::make_tuple(ret, x, y, s, t, dst, icurve);
        } 
        return data;
      }, 
        py::arg("qx"), py::arg("qy"), py::arg("icurve_begin"), py::arg("icurve_end"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the set of curve
        inside the range of segment ``icurve_begin..icurve_end``.
        Vectorial version.

        There are 7 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve
         7. **index** of the closest segment 

        :param List[float] qx: x coordinates of the point
        :param List[float] qy: y coordinates of the point
        :param int icurve_begin: first segment of the range
        :param int icurve_end: last segment of the range
        :return: a tuple of results as described
        :rtype: List[Tuple[int, float, float, float, float, float, int]]
      )S")

      .def("closestPointInRange_ISO", [](const ClothoidList & self, real_type qx, real_type qy, int_type icurve_begin, int_type icurve_end) {
        real_type x, y, s, t, dst;
        int_type icurve;
        int_type ret = self.closestPointInRange_ISO(qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve);
        return std::make_tuple(ret, x, y, s, t, dst, icurve);
      }, 
        py::arg("qx"), py::arg("qy"), py::arg("icurve_begin"), py::arg("icurve_end"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the set of curve
        inside the range of segment ``icurve_begin..icurve_end``.
        It uses the ISO reference frame

        There are 7 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve
         7. **index** of the closest segment 

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :param int icurve_begin: first segment of the range
        :param int icurve_end: last segment of the range
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, float, int]
      )S")

      .def("closestPointInRange_SAE", [](const ClothoidList & self, real_type qx, real_type qy, int_type icurve_begin, int_type icurve_end) {
        real_type x, y, s, t, dst;
        int_type icurve;
        int_type ret = self.closestPointInRange_SAE(qx, qy, icurve_begin, icurve_end, x, y, s, t, dst, icurve);
        return std::make_tuple(ret, x, y, s, t, dst, icurve);
      }, 
        py::arg("qx"), py::arg("qy"), py::arg("icurve_begin"), py::arg("icurve_end"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the set of curve
        inside the range of segment ``icurve_begin..icurve_end``.
        It uses the SAE reference frame.

        There are 7 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve
         7. **index** of the closest segment 

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :param int icurve_begin: first segment of the range
        :param int icurve_end: last segment of the range
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, float, int]
      )S")

      .def("closestPointInSRange", [](const ClothoidList & self, real_type qx, real_type qy, real_type s_begin, real_type s_end) {
        real_type x, y, s, t, dst;
        int_type icurve;
        int_type ret = self.closestPointInSRange_ISO(qx, qy, s_begin, s_begin, x, y, s, t, dst, icurve);
        return std::make_tuple(ret, x, y, s, t, dst, icurve);
      }, 
        py::arg("qx"), py::arg("qy"), py::arg("s_begin"), py::arg("s_end"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the set of curve
        inside the range provided as curvilinear abscissa :math:`[s_{start}, s_{end}]`

        There are 7 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve
         7. **index** of the closest segment 

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :param float s_begin: initial curvilinear abscissa
        :param float s_end: final curvilinear abscissa
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, float, int]
      )S")

      .def("closestPointInSRange", [](const ClothoidList & self, const std::vector<real_type> & qx, 
        const std::vector<real_type> & qy, real_type s_begin, real_type s_end) {
        const size_t n = std::min(qx.size(), qy.size());
        std::vector<std::tuple<int_type, real_type, real_type, real_type, real_type, real_type, int_type>> data(n);
        for (size_t i = 0; i < n; i++) {
          real_type x, y, s, t, dst;
          int_type icurve;
          int_type ret = self.closestPointInSRange_ISO(qx[i], qy[i], s_begin, s_end, x, y, s, t, dst, icurve);
          data[i] = std::make_tuple(ret, x, y, s, t, dst, icurve);
        } 
        return data;
      }, py::arg("qx"), py::arg("qy"), py::arg("s_begin"), py::arg("s_end"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the set of curve
        inside the range provided as curvilinear abscissa :math:`[s_{start}, s_{end}]`

        There are 7 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve
         7. **index** of the closest segment 

        :param List[float] qx: x coordinates of the point
        :param List[float] qy: y coordinates of the point
        :param float s_begin: initial curvilinear abscissa
        :param float s_end: final curvilinear abscissa
        :return: a tuple of results as described
        :rtype: List[Tuple[int, float, float, float, float, float, int]]
      )S")

      .def("closestPointInSRange_ISO", [](const ClothoidList & self, real_type qx, real_type qy, real_type s_begin, real_type s_end) {
        real_type x, y, s, t, dst;
        int_type icurve;
        int_type ret = self.closestPointInSRange_ISO(qx, qy, s_begin, s_begin, x, y, s, t, dst, icurve);
        return std::make_tuple(ret, x, y, s, t, dst, icurve);
      }, 
        py::arg("qx"), py::arg("qy"), py::arg("s_begin"), py::arg("s_end"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the set of curve
        inside the range provided as curvilinear abscissa :math:`[s_{start}, s_{end}]`
        It uses the ISO reference frame

        There are 7 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve
         7. **index** of the closest segment 

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :param float s_begin: initial curvilinear abscissa
        :param float s_end: final curvilinear abscissa
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, float, int]
      )S")

      .def("closestPointInSRange_SAE", [](const ClothoidList & self, real_type qx, real_type qy, real_type s_begin, real_type s_end) {
        real_type x, y, s, t, dst;
        int_type icurve;
        int_type ret = self.closestPointInSRange_SAE(qx, qy, s_begin, s_begin, x, y, s, t, dst, icurve);
        return std::make_tuple(ret, x, y, s, t, dst, icurve);
      }, 
        py::arg("qx"), py::arg("qy"), py::arg("s_begin"), py::arg("s_end"),
      R"S(
        Given a point :math:`(q_x, q_y)`, finds the closest point on the set of curve
        inside the range provided as curvilinear abscissa :math:`[s_{start}, s_{end}]`.
        It uses the SAE reference frame

        There are 7 values in the response tuple:

         1. result of the projection. If the value is 1 the projection is unique and 
            orthogonal, if the value is 0 there is more than one orthogonal projection
            and only the first is returned, if the value is -1 the minimum distance point 
            has a projection which is not orthogonal
         2. **x** coordinate of the point on the clothoid
         3. **y** coordinate of the point on the clothoid
         4. **s** curvilinear abscissa of the point
         5. **t** normal distance of the point on curvilinear abscissa
         6. **dst** distance of the point from the curve
         7. **index** of the closest segment 

        :param float qx: x coordinates of the point
        :param float qy: y coordinates of the point
        :param float s_begin: initial curvilinear abscissa
        :param float s_end: final curvilinear abscissa
        :return: a tuple of results as described
        :rtype: Tuple[int, float, float, float, float, float, int]
      )S")

      .def("build_AABBtree_ISO", &ClothoidList::build_AABBtree_ISO,
        py::arg("offs"), py::arg("max_angle") = Utils::m_pi/6, py::arg("max_size") = 1e100,
      R"S(
        Build the internal AABB tree of the clothoid list with offset (ISO)
        
        :param float offs:      curve offset
        :param float max_angle: maximum angle variation of the curve covered by a triangle
        :param float max_size:  maximum admissible size of the covering tirnagles
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("getSK", [](const ClothoidList & self) {
        std::vector<real_type> s, k;
        self.getSK(s, k);
        return std::make_tuple(s, k);
      },
      R"S(
        Returns the list of nodes angles and curvatures of the clothoid list

        :return: s and curvatures in nodes as lists
        :rtype: Tuple[List[float], List[float]]
      )S")

      .def("getSTK", [](const ClothoidList & self) {
        std::vector<real_type> s, t, k;
        self.getSTK(s, t, k);
        return std::make_tuple(s, t, k);
      },
      R"S(
        Returns the list of nodes angles and curvatures of the biarc

        :return: s, angles and curvatures in nodes as lists
        :rtype: Tuple[List[float], List[float], List[float]]
      )S")

      .def("getXY", [](const ClothoidList & self) {
        const size_t n = self.numSegments();
        // Avoids segmentation fault for empty list
        if (!n) {
          return make_tuple(std::vector<real_type>(), 
                            std::vector<real_type>());
        }
        std::vector<real_type> x(n+1), y(n+1);
        self.getXY(x.data(), y.data());
        return std::make_tuple(x, y);
      },
      R"S(
        Returns the list of nodes coordinates

        :return: nodes coordinates, a list of x and y values
        :rtype: Tuple[List[float], List[float]]
      )S")

      .def("getDeltaTheta", [](const ClothoidList & self) {
        const size_t n = self.numSegments();
        if (!n) { return std::vector<real_type>(); }
        std::vector<real_type> deltaTheta(n - 1);
        if (n - 1) {
          self.getDeltaTheta(deltaTheta.data());
        }
        return deltaTheta;
      },
      R"S(
        Returns the delta of thetas in nodes

        :return: delta theta list
        :rtype: List[float]
      )S")

      .def("getDeltaKappa", [](const ClothoidList & self) {
        const size_t n = self.numSegments();
        if (!n) { return std::vector<real_type>(); }
        std::vector<real_type> deltaKappa(n - 1);
        if (n - 1) {
          self.getDeltaKappa(deltaKappa.data());
        }
        return deltaKappa;
      },
      R"S(
        Returns the delta of curvature in nodes

        :return: delta curvature list
        :rtype: List[float]
      )S")

      .def("findST1", [](const ClothoidList & self, real_type x, real_type y) {
        real_type s, t;
        int_type idx = self.findST1(x, y, s, t);
        return std::make_tuple(idx, s, t);
      }, py::arg("x"), py::arg("y"),
      R"S(
        Find parametric coordinate :math:`(s, t)` given a point. With respect to the 
        classic `findST`, this version returns as first return value the index
        of the segment.

        :param float x: **x** coordinates of the point
        :param float y: **y** coordinates of the point
        :return: a tuple with index, **s** coordinates and **t** coordinates
        :rtype: Tuple[int, float, float]
      )S")

      .def("findST1", [](const ClothoidList & self, const std::vector<real_type> & x, const std::vector<real_type> & y) {
        const size_t n = std::min(x.size(), y.size());
        std::vector<real_type> s(n), t(n);
        std::vector<int_type> idx(n);
        for (size_t i = 0; i < n; i++) {
          idx[i] = self.findST1(x[i], y[i], s[i], t[i]);
        }
        return std::make_tuple(idx, s, t);
      }, py::arg("x"), py::arg("y"),
      R"S(
        Find parametric coordinate :math:`(s, t)` given a point. With respect to the 
        classic `findST`, this version returns as first return value the index
        of the segment.
        Vectorial version

        :param List[float] x: **x** coordinates of the point
        :param List[float] y: **y** coordinates of the point
        :return: a tuple with index, **s** coordinates and **t** coordinates
        :rtype: Tuple[List[int], List[float], List[float]]
      )S")

      .def("findST1", [](const ClothoidList & self, int_type ibegin, int_type iend, real_type x, real_type y) {
        real_type s, t;
        int_type idx = self.findST1(ibegin, iend, x, y, s, t);
        return std::make_tuple(idx, s, t);
      }, py::arg("ibegin"), py::arg("iend"), py::arg("x"), py::arg("y"),
      R"S(
        Find parametric coordinate :math:`(s, t)` given a point. With respect to the 
        classic `findST`, this version returns as first return value the index
        of the segment. The search is done only between elements ``ibegin..iend``.

        :param int ibegin: initial index
        :param int iend: final index
        :param float x: **x** coordinates of the point
        :param float y: **y** coordinates of the point
        :return: a tuple with index, **s** coordinates and **t** coordinates
        :rtype: Tuple[int, float, float]
      )S")

      .def("findST1", [](const ClothoidList & self, int_type ibegin, int_type iend, const std::vector<real_type> & x, const std::vector<real_type> & y) {
        const size_t n = std::min(x.size(), y.size());
        std::vector<real_type> s(n), t(n);
        std::vector<int_type> idx(n);
        for (size_t i = 0; i < n; i++) {
          idx[i] = self.findST1(ibegin, iend, x[i], y[i], s[i], t[i]);
        }
        return std::make_tuple(idx, s, t);
      }, py::arg("ibegin"), py::arg("iend"), py::arg("x"), py::arg("y"),
      R"S(
        Find parametric coordinate :math:`(s, t)` given a point. With respect to the 
        classic `findST`, this version returns as first return value the index
        of the segment. The search is done only between elements ``ibegin..iend``.
        Vectorial version.

        :param int ibegin: initial index
        :param int iend: final index
        :param List[float] x: **x** coordinates of the point
        :param List[float] y: **y** coordinates of the point
        :return: a tuple with index, **s** coordinates and **t** coordinates
        :rtype: Tuple[List[int], List[float], List[float]]
      )S")

      .def("save", [](const ClothoidList & self) {
        std::ostringstream str;
        self.save(str);
        return str.str();
      },
      R"S(
        Export the current clothoid list in a CSV

        :return: the CSV string
        :rtype: str
      )S")

      .def("export_table", [](const ClothoidList & self) {
        std::ostringstream str;
        self.export_table(str);
        return str.str();
      },
      R"S(
        Export the current clothoid list in a table

        :return: the table string
        :rtype: str
      )S")

      .def("load", [](ClothoidList & self, const std::string & csv, real_type epsi) {
        std::istringstream str(csv);
        self.load(str, epsi);
      }, py::arg("csv"), py::arg("espi") = 1e-8,
      R"S(
        Load the clothoid from a csv. Data is assumed to be saved as follows:

            # x y theta kappa
            x0 y0 theta0 kappa0
            x1 y1 theta1 kappa1
            ...
            xn yn thetan kappan

        like data exported in save.

        :param str csv: the stream of data
        :param float epsi: precision
        :return: nothing, works in place
        :rtype: NoneType
      )S")

      .def("as_list", [](const ClothoidList & self) {
        const size_t n = self.numSegments();
        std::vector<ClothoidCurve> list;
        list.reserve(n);
        for (size_t i = 0; i < n; i++) {
          list.push_back(self.get(i));
        }
        return list;
      }, 
      R"S(
        Return the clothoid list as a python list

        :return: a list with ClothoidCurve objects
        :rtype: List[ClothoidCurve]
      )S");
    }
  }
}