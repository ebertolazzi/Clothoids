/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2025                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Università degli Studi di Trento                                    |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/*
 *  Note:
 *
 *    Most of the problems were originally part of ACM algorithm 566,
 *    and were used to test the MINPACK package.
 *
 *  Author:
 *
 *    John Burkardt
 *
 *  Reference:
 *

@book{Dennis:1996,
  author    = {Dennis, J. and Schnabel, R.},
  title     = {Numerical Methods for Unconstrained Optimization and Nonlinear
Equations}, publisher = {Society for Industrial and Applied Mathematics}, year
= {1996}, doi       = {10.1137/1.9781611971200},
}

 *
 *    N. de Villiers and D. Glasser,
 *    A continuation method for nonlinear regression,
 *    SIAM Journal of Numerical Analysis,
 *    Volume 18, 1981, pages 1139-1154.
 *
 *    C. Fraley,
 *    Solution of nonlinear least-squares problems,
 *    Technical Report STAN-CS-1165,
 *    Computer Science Department,
 *    Stanford University, 1987.
 *
 *    C. Fraley,
 *    Software performance on nonlinear least-squares problems,
 *    Technical Report SOL 88-17,
 *    Systems Optimization Laboratory,
 *    Department of Operations Research,
 *    Stanford University, 1988.
 *
 *    J. J. McKeown,
 *    Specialized versus general-purpose algorithms for functions
 *    that are sums of squared terms,
 *    Math. Prog.,
 *    Volume 9, 1975, pages 57-68.
 *
 *    J. J. McKeown,
 *    On algorithms for sums of squares problems,
 *    in Towards Global Optimization,
 *    L. C. W. Dixon and G. Szego (editors),
 *    North-Holland (1975b), pages 229-257.
 *
 *
 *    D. E. Salane,
 *    A continuation approach for solving large residual
 *    nonlinear least squares problems,
 *    SIAM Journal of Scientific and Statistical Computing,
 *    Volume 8, 1987, pages 655-671.
 */

#include "Utils_nonlinear_system.hh"

namespace Utils
{

  using real_type = double;

  static inline real_type power2( real_type a )
  {
    return a * a;
  }
  static inline real_type power3( real_type a )
  {
    return a * a * a;
  }
  static inline real_type power4( real_type a )
  {
    return a * a * a * a;
  }
  static inline real_type power5( real_type a )
  {
    return a * a * a * a * a;
  }
  static inline real_type power6( real_type a )
  {
    return a * a * a * a * a * a;
  }

#include "Utils_nonlinear_system/ArtificialTestOfNowakAndWeimann.hxx"
#include "Utils_nonlinear_system/BUNLSI.hxx"
#include "Utils_nonlinear_system/BadlyScaledAugmentedPowellFunction.hxx"
#include "Utils_nonlinear_system/Beale.hxx"
#include "Utils_nonlinear_system/Bertolazzi.hxx"
#include "Utils_nonlinear_system/BiggsEXPfunctions.hxx"
#include "Utils_nonlinear_system/BoggsFunction.hxx"
#include "Utils_nonlinear_system/Bohachevsky.hxx"
#include "Utils_nonlinear_system/Box3.hxx"
#include "Utils_nonlinear_system/BoxAndBettsExponentialQuadraticSum.hxx"
#include "Utils_nonlinear_system/BoxProblem.hxx"
#include "Utils_nonlinear_system/BraninRCOS.hxx"
#include "Utils_nonlinear_system/BrownAlmostLinearFunction.hxx"
#include "Utils_nonlinear_system/BrownAndConteFunction.hxx"
#include "Utils_nonlinear_system/BrownAndDennis.hxx"
#include "Utils_nonlinear_system/BrownAndGearhartFunction.hxx"
#include "Utils_nonlinear_system/BrownBadlyScaled.hxx"
#include "Utils_nonlinear_system/BrownFunction.hxx"
#include "Utils_nonlinear_system/BroydenBandedFunction.hxx"
#include "Utils_nonlinear_system/BroydenTridiagonalFunction.hxx"
#include "Utils_nonlinear_system/BurdenAndFaires.hxx"
#include "Utils_nonlinear_system/Chandrasekhar.hxx"
#include "Utils_nonlinear_system/ChebyquadFunction.hxx"
#include "Utils_nonlinear_system/ChemicalEquilibriumApplication.hxx"
#include "Utils_nonlinear_system/CliffFunction.hxx"
#include "Utils_nonlinear_system/Colville.hxx"
#include "Utils_nonlinear_system/CombustionApplication.hxx"
#include "Utils_nonlinear_system/ComplementaryFunction.hxx"
#include "Utils_nonlinear_system/CompressibilityFactorFromTheRKequation.hxx"
#include "Utils_nonlinear_system/CountercurrentReactorsProblem.hxx"
#include "Utils_nonlinear_system/CraggAndLevyProblem.hxx"
#include "Utils_nonlinear_system/CubeFunction.hxx"
#include "Utils_nonlinear_system/DarvishiBarati.hxx"
#include "Utils_nonlinear_system/DeVilliersGlasser.hxx"
#include "Utils_nonlinear_system/DeistSefor.hxx"
#include "Utils_nonlinear_system/DennisAndGay.hxx"
#include "Utils_nonlinear_system/DennisAndSchnabel2x2example.hxx"
#include "Utils_nonlinear_system/DiagonalFunctionMulQO.hxx"
#include "Utils_nonlinear_system/DiscreteBoundaryValueFunction.hxx"
#include "Utils_nonlinear_system/DiscreteIntegralEquationFunction.hxx"
#include "Utils_nonlinear_system/DixonFunction.hxx"
#include "Utils_nonlinear_system/Easom.hxx"
#include "Utils_nonlinear_system/EsterificReaction.hxx"
#include "Utils_nonlinear_system/ExponentialFunction.hxx"
#include "Utils_nonlinear_system/ExponentialSine.hxx"
#include "Utils_nonlinear_system/ExtendedEigerSikorskiStenger.hxx"
#include "Utils_nonlinear_system/ExtendedKearfottFunction.hxx"
#include "Utils_nonlinear_system/ExtendedPowellSingularFunction.hxx"
#include "Utils_nonlinear_system/FreudensteinRothFunction.hxx"
#include "Utils_nonlinear_system/Function15.hxx"
#include "Utils_nonlinear_system/Function18.hxx"
#include "Utils_nonlinear_system/Function21.hxx"
#include "Utils_nonlinear_system/Function27.hxx"
#include "Utils_nonlinear_system/Gauss.hxx"
#include "Utils_nonlinear_system/GeneralizedRosenbrock.hxx"
#include "Utils_nonlinear_system/GeometricProgrammingFunction.hxx"
#include "Utils_nonlinear_system/GheriMancino.hxx"
#include "Utils_nonlinear_system/GoldsteinPrice.hxx"
#include "Utils_nonlinear_system/GregoryAndKarney.hxx"
#include "Utils_nonlinear_system/GriewankFunction.hxx"
#include "Utils_nonlinear_system/Gulf.hxx"
#include "Utils_nonlinear_system/HAS.hxx"
#include "Utils_nonlinear_system/Hammarling.hxx"
#include "Utils_nonlinear_system/HanSunHan.hxx"
#include "Utils_nonlinear_system/HanbookFunction.hxx"
#include "Utils_nonlinear_system/HelicalValleyFunction.hxx"
#include "Utils_nonlinear_system/HiebertChem.hxx"
#include "Utils_nonlinear_system/Hilbert.hxx"
#include "Utils_nonlinear_system/Himmelblau.hxx"
#include "Utils_nonlinear_system/InfRefluxFunction.hxx"
#include "Utils_nonlinear_system/IntervalArithmeticBenchmarks.hxx"
#include "Utils_nonlinear_system/JennrichAndSampsonFunction.hxx"
#include "Utils_nonlinear_system/KelleyFunction.hxx"
#include "Utils_nonlinear_system/KinematicApplication.hxx"
#include "Utils_nonlinear_system/Leon.hxx"
#include "Utils_nonlinear_system/LinearFunctionFullRank.hxx"
#include "Utils_nonlinear_system/LinearFunctionRank1.hxx"
#include "Utils_nonlinear_system/LogarithmicFunction.hxx"
#include "Utils_nonlinear_system/McCormicFunction.hxx"
#include "Utils_nonlinear_system/McKinnon.hxx"
#include "Utils_nonlinear_system/MexicanHatFunction.hxx"
#include "Utils_nonlinear_system/MieleAndCantrellFunction.hxx"
#include "Utils_nonlinear_system/NonlinearIntegralEquations.hxx"
#include "Utils_nonlinear_system/Order10to11function.hxx"
#include "Utils_nonlinear_system/PavianiFunction.hxx"
#include "Utils_nonlinear_system/Penalty.hxx"
#include "Utils_nonlinear_system/Powell3D.hxx"
#include "Utils_nonlinear_system/PowellBadlyScaledFunction.hxx"
#include "Utils_nonlinear_system/PowellQuarticFunction.hxx"
#include "Utils_nonlinear_system/RooseKullaLombMeressoo.hxx"
#include "Utils_nonlinear_system/SIRtest.hxx"
#include "Utils_nonlinear_system/SSTnonlinearityTerm.hxx"
#include "Utils_nonlinear_system/SampleProblem18.hxx"
#include "Utils_nonlinear_system/SampleProblem19.hxx"
#include "Utils_nonlinear_system/ScalarProblem.hxx"
#include "Utils_nonlinear_system/Schaffer.hxx"
#include "Utils_nonlinear_system/SchubertBroydenFunction.hxx"
#include "Utils_nonlinear_system/Semiconductor2D.hxx"
#include "Utils_nonlinear_system/Shacham.hxx"
#include "Utils_nonlinear_system/Shekel.hxx"
#include "Utils_nonlinear_system/ShenYpma.hxx"
#include "Utils_nonlinear_system/Shubert.hxx"
#include "Utils_nonlinear_system/SingularFunction.hxx"
#include "Utils_nonlinear_system/SingularSystem.hxx"
#include "Utils_nonlinear_system/SixHumpCamelBackFunction.hxx"
#include "Utils_nonlinear_system/SoniaKrzyworzcka.hxx"
#include "Utils_nonlinear_system/Spedicato.hxx"
#include "Utils_nonlinear_system/StrictlyConvexFunction.hxx"
#include "Utils_nonlinear_system/Toint.hxx"
#include "Utils_nonlinear_system/TridimensionalValley.hxx"
#include "Utils_nonlinear_system/TrigonometricExponentialSystem.hxx"
#include "Utils_nonlinear_system/TrigonometricFunction.hxx"
#include "Utils_nonlinear_system/TroeschFunction.hxx"
#include "Utils_nonlinear_system/TwoPointBoundaryValueProblem.hxx"
#include "Utils_nonlinear_system/VariablyDimensionedFunction.hxx"
#include "Utils_nonlinear_system/WatsonFunction.hxx"
#include "Utils_nonlinear_system/Weibull.hxx"
#include "Utils_nonlinear_system/WoodFunction.hxx"
#include "Utils_nonlinear_system/XiaoYin.hxx"
#include "Utils_nonlinear_system/YixunShi.hxx"
#include "Utils_nonlinear_system/ZeroJacobianFunction.hxx"

  std::vector<NonlinearSystem *> nonlinear_system_tests;
  std::map<string, unsigned>     nonlinear_system_tests_map;

  void init_nonlinear_system_tests()
  {
    nonlinear_system_tests.push_back( new ArtificialTestOfNowakAndWeimann() );
    nonlinear_system_tests.push_back( new BadlyScaledAugmentedPowellFunction( 3 ) );
    nonlinear_system_tests.push_back( new BadlyScaledAugmentedPowellFunction( 30 ) );
    nonlinear_system_tests.push_back( new BadlyScaledAugmentedPowellFunction( 300 ) );
    // nonlinear_system_tests.push_back(new BadlyScaledAugmentedPowellFunction(3000));
    nonlinear_system_tests.push_back( new Beale() );
    nonlinear_system_tests.push_back( new BertolazziRootPlusSquare() );
    nonlinear_system_tests.push_back( new BertolazziAtanPlusQuadratic() );
    nonlinear_system_tests.push_back( new BertolazziHard() );
    nonlinear_system_tests.push_back( new BertolazziSingleEQ() );

    nonlinear_system_tests.push_back( new BiggsEXP2function() );
    nonlinear_system_tests.push_back( new BiggsEXP3function() );
    nonlinear_system_tests.push_back( new BiggsEXP4function() );
    nonlinear_system_tests.push_back( new BiggsEXP5function() );
    nonlinear_system_tests.push_back( new BiggsEXP6function() );
    nonlinear_system_tests.push_back( new BoggsFunction() );
    nonlinear_system_tests.push_back( new BohachevskyN1() );
    nonlinear_system_tests.push_back( new BohachevskyN2() );
    nonlinear_system_tests.push_back( new BohachevskyN3() );

    nonlinear_system_tests.push_back( new BoxAndBettsExponentialQuadraticSum() );
    nonlinear_system_tests.push_back( new BoxProblem() );
    nonlinear_system_tests.push_back( new Box3() );
    nonlinear_system_tests.push_back( new BraninRCOS() );
    nonlinear_system_tests.push_back( new BrownAlmostLinearFunction( 5 ) );
    nonlinear_system_tests.push_back( new BrownAlmostLinearFunction( 15 ) );
    nonlinear_system_tests.push_back( new BrownAlmostLinearFunction( 25 ) );
    nonlinear_system_tests.push_back( new BrownAndConteFunction() );
    nonlinear_system_tests.push_back( new BrownAndDennis() );
    nonlinear_system_tests.push_back( new BrownAndGearhartFunction() );
    nonlinear_system_tests.push_back( new BrownBadlyScaled() );
    nonlinear_system_tests.push_back( new BrownFunction() );
    nonlinear_system_tests.push_back( new BroydenBandedFunction() );

    nonlinear_system_tests.push_back( new BroydenTridiagonalFunction( 0.1, 1, 5 ) );
    nonlinear_system_tests.push_back( new BroydenTridiagonalFunction( 0.1, 1, 10 ) );
    nonlinear_system_tests.push_back( new BroydenTridiagonalFunction( 0.1, 1, 500 ) );
    nonlinear_system_tests.push_back( new BroydenTridiagonalFunction( 0.5, 1, 5 ) );
    nonlinear_system_tests.push_back( new BroydenTridiagonalFunction( 0.5, 1, 10 ) );
    nonlinear_system_tests.push_back( new BroydenTridiagonalFunction( 0.5, 1, 500 ) );

    nonlinear_system_tests.push_back( new BUNLSI5() );
    nonlinear_system_tests.push_back( new BUNLSI6() );

    nonlinear_system_tests.push_back( new BurdenAndFaires() );

    nonlinear_system_tests.push_back( new Chandrasekhar( 0.9999, 10 ) );
    nonlinear_system_tests.push_back( new Chandrasekhar( 0.9999, 50 ) );
    nonlinear_system_tests.push_back( new Chandrasekhar( 0.9, 10 ) );
    nonlinear_system_tests.push_back( new Chandrasekhar( 0.9, 50 ) );

    nonlinear_system_tests.push_back( new ChebyquadFunction( 1 ) );
    nonlinear_system_tests.push_back( new ChebyquadFunction( 2 ) );
    nonlinear_system_tests.push_back( new ChebyquadFunction( 3 ) );
    nonlinear_system_tests.push_back( new ChebyquadFunction( 4 ) );
    nonlinear_system_tests.push_back( new ChebyquadFunction( 5 ) );
    nonlinear_system_tests.push_back( new ChebyquadFunction( 6 ) );
    nonlinear_system_tests.push_back( new ChebyquadFunction( 7 ) );
    // nonlinear_system_tests.push_back( new ChebyquadFunction(8) );
    nonlinear_system_tests.push_back( new ChebyquadFunction( 9 ) );

    nonlinear_system_tests.push_back( new ChemicalEquilibriumApplication() );
    nonlinear_system_tests.push_back( new ChemicalEquilibriumPartialMethaneOxidation() );
    nonlinear_system_tests.push_back( new ChemicalReactorEquilibriumConversion() );
    nonlinear_system_tests.push_back( new ChemicalReactorSteadyState() );

    nonlinear_system_tests.push_back( new CliffFunction() );
    nonlinear_system_tests.push_back( new Colville() );
    nonlinear_system_tests.push_back( new CombustionApplication() );

    nonlinear_system_tests.push_back( new ComplementaryFunction( 2 ) );
    nonlinear_system_tests.push_back( new ComplementaryFunction( 16 ) );
    nonlinear_system_tests.push_back( new ComplementaryFunction( 128 ) );

    nonlinear_system_tests.push_back( new CompressibilityFactorFromTheRKequation() );

    nonlinear_system_tests.push_back( new CountercurrentReactorsProblem1( 6 ) );
    nonlinear_system_tests.push_back( new CountercurrentReactorsProblem1( 12 ) );
    nonlinear_system_tests.push_back( new CountercurrentReactorsProblem1( 50 ) );
    nonlinear_system_tests.push_back( new CountercurrentReactorsProblem2( 6 ) );
    nonlinear_system_tests.push_back( new CountercurrentReactorsProblem2( 12 ) );
    nonlinear_system_tests.push_back( new CountercurrentReactorsProblem2( 50 ) );
    nonlinear_system_tests.push_back( new CraggAndLevyProblem() );
    nonlinear_system_tests.push_back( new CubeFunction() );

    nonlinear_system_tests.push_back( new CutlipsSteadyStateForReactionRateEquations( 0 ) );
    nonlinear_system_tests.push_back( new CutlipsSteadyStateForReactionRateEquations( 1 ) );
    nonlinear_system_tests.push_back( new CutlipsSteadyStateForReactionRateEquations( 2 ) );

    nonlinear_system_tests.push_back( new DarvishiBarati() );

    nonlinear_system_tests.push_back( new DennisAndGay6eqN1() );
    nonlinear_system_tests.push_back( new DennisAndGay6eqN2() );
    nonlinear_system_tests.push_back( new DennisAndGay6eqN3() );
    nonlinear_system_tests.push_back( new DennisAndGay6eqN4() );
    nonlinear_system_tests.push_back( new DennisAndGay6eqN5() );
    nonlinear_system_tests.push_back( new DennisAndGay8eqN1() );
    nonlinear_system_tests.push_back( new DennisAndGay8eqN2() );
    nonlinear_system_tests.push_back( new DennisAndGay8eqN3() );
    nonlinear_system_tests.push_back( new DennisAndGay8eqN4() );
    nonlinear_system_tests.push_back( new DennisAndGay8eqN5() );

    nonlinear_system_tests.push_back( new DennisAndSchnabel2x2example() );

    nonlinear_system_tests.push_back( new DeVilliersGlasser01() );
    nonlinear_system_tests.push_back( new DeVilliersGlasser02() );

    nonlinear_system_tests.push_back( new DiagonalFunctionMulQO( 9 ) );
    nonlinear_system_tests.push_back( new DiagonalFunctionMulQO( 27 ) );

    nonlinear_system_tests.push_back( new DiscreteBoundaryValueFunction( 10 ) );
    nonlinear_system_tests.push_back( new DiscreteBoundaryValueFunction( 50 ) );
    nonlinear_system_tests.push_back( new DiscreteBoundaryValueFunction( 100 ) );
    nonlinear_system_tests.push_back( new DiscreteBoundaryValueFunction( 500 ) );
    // nonlinear_system_tests.push_back(new DiscreteBoundaryValueFunction(5000));

    nonlinear_system_tests.push_back( new DiscreteIntegralEquationFunction( 2 ) );
    nonlinear_system_tests.push_back( new DiscreteIntegralEquationFunction( 5 ) );
    nonlinear_system_tests.push_back( new DiscreteIntegralEquationFunction( 10 ) );
    nonlinear_system_tests.push_back( new DiscreteIntegralEquationFunction( 100 ) );

    nonlinear_system_tests.push_back( new DixonFunction( 80 ) );
    nonlinear_system_tests.push_back( new DixonFunction( 2000 ) );
    // nonlinear_system_tests.push_back(new DixonFunction(5000));

    nonlinear_system_tests.push_back( new Easom() );
    nonlinear_system_tests.push_back( new EsterificReaction() );

    nonlinear_system_tests.push_back( new ExponentialFunction1( 2 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction1( 10 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction1( 50 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction1( 500 ) );
    // nonlinear_system_tests.push_back(new ExponentialFunction1(5000));
    nonlinear_system_tests.push_back( new ExponentialFunction2( 2 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction2( 10 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction2( 50 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction2( 500 ) );
    // nonlinear_system_tests.push_back(new ExponentialFunction2(5000));
    nonlinear_system_tests.push_back( new ExponentialFunction3( 2 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction3( 10 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction3( 50 ) );
    nonlinear_system_tests.push_back( new ExponentialFunction3( 500 ) );
    // nonlinear_system_tests.push_back(new ExponentialFunction3(5000));

    nonlinear_system_tests.push_back( new ExponentialSine() );
    nonlinear_system_tests.push_back( new ExtendedEigerSikorskiStenger() );
    nonlinear_system_tests.push_back( new ExtendedKearfottFunction() );
    nonlinear_system_tests.push_back( new ExtendedPowellSingularFunction() );

    nonlinear_system_tests.push_back( new FractionalConversionInAchemicalReactor() );
    nonlinear_system_tests.push_back( new FractionalConversionInAchemicalReactor2() );

    nonlinear_system_tests.push_back( new FreudensteinRothFunction() );

    nonlinear_system_tests.push_back( new Function15( 10 ) );
    nonlinear_system_tests.push_back( new Function15( 50 ) );
    nonlinear_system_tests.push_back( new Function15( 100 ) );

    nonlinear_system_tests.push_back( new Function18( 3 ) );
    nonlinear_system_tests.push_back( new Function18( 9 ) );
    nonlinear_system_tests.push_back( new Function18( 27 ) );

    nonlinear_system_tests.push_back( new Function21( 3 ) );
    nonlinear_system_tests.push_back( new Function21( 6 ) );
    nonlinear_system_tests.push_back( new Function21( 36 ) );

    nonlinear_system_tests.push_back( new Function27( 2 ) );
    nonlinear_system_tests.push_back( new Function27( 10 ) );
    nonlinear_system_tests.push_back( new Function27( 100 ) );

    nonlinear_system_tests.push_back( new Gauss() );

    nonlinear_system_tests.push_back( new GeneralizedRosenbrock( 2 ) );
    nonlinear_system_tests.push_back( new GeneralizedRosenbrock( 10 ) );
    nonlinear_system_tests.push_back( new GeneralizedRosenbrock( 50 ) );
    nonlinear_system_tests.push_back( new GeneralizedRosenbrock( 500 ) );

    nonlinear_system_tests.push_back( new GeometricProgrammingFunction( 5 ) );
    nonlinear_system_tests.push_back( new GeometricProgrammingFunction( 10 ) );
    nonlinear_system_tests.push_back( new GeometricProgrammingFunction( 50 ) );
    nonlinear_system_tests.push_back( new GeometricProgrammingFunction( 100 ) );

    nonlinear_system_tests.push_back( new GheriMancino( 10 ) );
    nonlinear_system_tests.push_back( new GheriMancino( 30 ) );
    nonlinear_system_tests.push_back( new GheriMancino( 100 ) );

    nonlinear_system_tests.push_back( new GoldsteinPrice() );
    nonlinear_system_tests.push_back( new GregoryAndKarney( 10 ) );
    nonlinear_system_tests.push_back( new GregoryAndKarney( 20 ) );
    nonlinear_system_tests.push_back( new GregoryAndKarney( 100 ) );

    nonlinear_system_tests.push_back( new GriewankFunction( 2 ) );
    nonlinear_system_tests.push_back( new GriewankFunction( 5 ) );
    nonlinear_system_tests.push_back( new GriewankFunction( 10 ) );

    nonlinear_system_tests.push_back( new Gulf() );

    nonlinear_system_tests.push_back( new Hammarling2x2matrixSquareRoot() );
    nonlinear_system_tests.push_back( new Hammarling3x3matrixSquareRootProblemN1() );
    nonlinear_system_tests.push_back( new Hammarling3x3matrixSquareRootProblemN2() );
    nonlinear_system_tests.push_back( new Hammarling3x3matrixSquareRootProblemN3() );

    nonlinear_system_tests.push_back( new HanbookFunction( 2 ) );
    nonlinear_system_tests.push_back( new HanbookFunction( 10 ) );
    nonlinear_system_tests.push_back( new HanbookFunction( 50 ) );

    nonlinear_system_tests.push_back( new HanSunHan() );

    nonlinear_system_tests.push_back( new HAS64( 1e2 ) );
    nonlinear_system_tests.push_back( new HAS64( 1e4 ) );
    nonlinear_system_tests.push_back( new HAS64( 1e6 ) );
    nonlinear_system_tests.push_back( new HAS64( 1e8 ) );
    nonlinear_system_tests.push_back( new HAS64( 1e10 ) );
    nonlinear_system_tests.push_back( new HAS93( 1e2 ) );
    nonlinear_system_tests.push_back( new HAS93( 1e4 ) );
    nonlinear_system_tests.push_back( new HAS93( 1e6 ) );
    nonlinear_system_tests.push_back( new HAS93( 1e8 ) );
    nonlinear_system_tests.push_back( new HAS93( 1e10 ) );
    nonlinear_system_tests.push_back( new HAS111() );

    nonlinear_system_tests.push_back( new HelicalValleyFunction() );

    nonlinear_system_tests.push_back( new HiebertChem10x10() );
    nonlinear_system_tests.push_back( new HiebertChem2x2() );
    nonlinear_system_tests.push_back( new HiebertChem6x6() );

    nonlinear_system_tests.push_back( new Hiebert3ChemicalEquilibriumProblem( 10 ) );
    nonlinear_system_tests.push_back( new Hiebert3ChemicalEquilibriumProblem( 40 ) );

    nonlinear_system_tests.push_back( new Hilbert( 4 ) );
    nonlinear_system_tests.push_back( new Hilbert( 8 ) );
    nonlinear_system_tests.push_back( new Hilbert( 32 ) );
    nonlinear_system_tests.push_back( new Hilbert( 64 ) );
    nonlinear_system_tests.push_back( new Hilbert( 256 ) );

    nonlinear_system_tests.push_back( new Himmelblau() );
    nonlinear_system_tests.push_back( new InfRefluxFunction() );
    nonlinear_system_tests.push_back( new IntervalArithmeticBenchmarks() );
    nonlinear_system_tests.push_back( new JennrichAndSampsonFunction() );
    nonlinear_system_tests.push_back( new KelleyFunction() );
    nonlinear_system_tests.push_back( new KinematicApplication() );
    nonlinear_system_tests.push_back( new Leon() );
    nonlinear_system_tests.push_back( new LinearFunctionFullRank() );
    nonlinear_system_tests.push_back( new LinearFunctionRank1() );

    nonlinear_system_tests.push_back( new LogarithmicFunction( 2 ) );
    nonlinear_system_tests.push_back( new LogarithmicFunction( 10 ) );
    nonlinear_system_tests.push_back( new LogarithmicFunction( 50 ) );
    nonlinear_system_tests.push_back( new LogarithmicFunction( 500 ) );
    // nonlinear_system_tests.push_back(new LogarithmicFunction(5000));

    nonlinear_system_tests.push_back( new McCormicFunction() );
    nonlinear_system_tests.push_back( new McKinnon() );

    nonlinear_system_tests.push_back( new MexicanHatFunction( 1e2 ) );
    nonlinear_system_tests.push_back( new MexicanHatFunction( 1e4 ) );
    nonlinear_system_tests.push_back( new MexicanHatFunction( 1e6 ) );
    nonlinear_system_tests.push_back( new MexicanHatFunction( 1e8 ) );

    nonlinear_system_tests.push_back( new MieleAndCantrellFunction() );

    nonlinear_system_tests.push_back( new ModelEquationsForTheCSTR() );
    nonlinear_system_tests.push_back( new ModelEquationsForCombustionOfPropane( 5 ) );
    nonlinear_system_tests.push_back( new ModelEquationsForCombustionOfPropane( 10 ) );

    nonlinear_system_tests.push_back( new NonlinearIntegralEquations() );
    nonlinear_system_tests.push_back( new Order10to11function() );

    nonlinear_system_tests.push_back( new PavianiFunction() );

    nonlinear_system_tests.push_back( new PenaltyIfunction( 10 ) );
    nonlinear_system_tests.push_back( new PenaltyIfunction( 50 ) );

    nonlinear_system_tests.push_back( new PenaltyN1( 2 ) );
    nonlinear_system_tests.push_back( new PenaltyN1( 10 ) );
    nonlinear_system_tests.push_back( new PenaltyN1( 50 ) );
    nonlinear_system_tests.push_back( new PenaltyN1( 200 ) );

    nonlinear_system_tests.push_back( new PenaltyN2( 2 ) );
    nonlinear_system_tests.push_back( new PenaltyN2( 10 ) );
    nonlinear_system_tests.push_back( new PenaltyN2( 50 ) );
    nonlinear_system_tests.push_back( new PenaltyN2( 200 ) );

    nonlinear_system_tests.push_back( new PipelineNetworkProblem() );
    nonlinear_system_tests.push_back( new PipelineNetworkProblem2() );

    nonlinear_system_tests.push_back( new PowellBadlyScaledFunction() );
    nonlinear_system_tests.push_back( new PowellQuarticFunction() );
    nonlinear_system_tests.push_back( new Powell3D() );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo129() );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo201( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo201( 100 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo201( 500 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo202( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo202( 100 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo203( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo203( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo203( 7 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo203( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo203( 15 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo203( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo203( 30 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo203( 50 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo204( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo204( 100 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo204( 500 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo205( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo205( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo205( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo205( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo205( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo206( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo206( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo206( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo206( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo206( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo207( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo207( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo207( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo207( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo207( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo208( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo208( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo208( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo208( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo208( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo209( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo209( 100 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo209( 500 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo210( 10 ) );  // no solution
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo211( 10 ) );  // no solution

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo212( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo212( 100 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo212( 500 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo213( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo213( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo213( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo213( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo213( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo214( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo214( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo214( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo214( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo214( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo215( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo215( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo215( 6 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo215( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo215( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo215( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo216( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo216( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo216( 7 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo216( 9 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo217( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo217( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo217( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo217( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo217( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo218( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo218( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo218( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo218( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo218( 30 ) );

    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo219( 2 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo219( 5 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo219( 10 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo219( 20 ) );
    nonlinear_system_tests.push_back( new RooseKullaLombMeressoo219( 30 ) );

    nonlinear_system_tests.push_back( new SampleProblem18() );
    nonlinear_system_tests.push_back( new SampleProblem19() );
    nonlinear_system_tests.push_back( new ScalarProblem() );

    nonlinear_system_tests.push_back( new SchafferF6() );
    nonlinear_system_tests.push_back( new SchafferF7() );

    nonlinear_system_tests.push_back( new SchubertBroydenFunction( 10 ) );
    nonlinear_system_tests.push_back( new SchubertBroydenFunction( 100 ) );
    nonlinear_system_tests.push_back( new SchubertBroydenFunction( 1000 ) );
    // nonlinear_system_tests.push_back(new SchubertBroydenFunction(5000));

    nonlinear_system_tests.push_back( new Semiconductor2D() );

    nonlinear_system_tests.push_back( new ShekelSQRN5() );
    nonlinear_system_tests.push_back( new ShekelSQRN7() );
    nonlinear_system_tests.push_back( new ShekelSQRN10() );

    nonlinear_system_tests.push_back( new ShenYpma5() );
    nonlinear_system_tests.push_back( new ShenYpma7() );
    nonlinear_system_tests.push_back( new ShenYpma8() );

    nonlinear_system_tests.push_back( new Shubert() );

    nonlinear_system_tests.push_back( new SingularFunction( 2 ) );
    nonlinear_system_tests.push_back( new SingularFunction( 10 ) );
    nonlinear_system_tests.push_back( new SingularFunction( 50 ) );
    nonlinear_system_tests.push_back( new SingularFunction( 500 ) );
    // nonlinear_system_tests.push_back(new SingularFunction(5000));

    nonlinear_system_tests.push_back( new SingularSystemA() );
    nonlinear_system_tests.push_back( new SingularSystemB() );
    nonlinear_system_tests.push_back( new SingularSystemC() );
    nonlinear_system_tests.push_back( new SingularSystemD() );
    nonlinear_system_tests.push_back( new SingularSystemE() );
    nonlinear_system_tests.push_back( new SingularSystemF() );

    nonlinear_system_tests.push_back( new SingularSystemP2() );
    nonlinear_system_tests.push_back( new SingularSystemP3() );
    nonlinear_system_tests.push_back( new SingularSystemP4() );
    nonlinear_system_tests.push_back( new SingularSystemP5() );
    nonlinear_system_tests.push_back( new SingularSystemP6() );
    nonlinear_system_tests.push_back( new SingularSystemP7() );
    nonlinear_system_tests.push_back( new SingularSystemP8() );
    nonlinear_system_tests.push_back( new SingularSystemP9() );

    nonlinear_system_tests.push_back( new SIRtest( 2 ) );
    nonlinear_system_tests.push_back( new SIRtest( 20 ) );
    nonlinear_system_tests.push_back( new SIRtest( 100 ) );

    nonlinear_system_tests.push_back( new SixHumpCamelBackFunction() );

    nonlinear_system_tests.push_back( new SoniaKrzyworzcka1() );
    nonlinear_system_tests.push_back( new SoniaKrzyworzcka2() );

    nonlinear_system_tests.push_back( new SpedicatoFunction17( 10 ) );
    nonlinear_system_tests.push_back( new SpedicatoFunction17( 50 ) );
    nonlinear_system_tests.push_back( new SpedicatoFunction17( 100 ) );
    nonlinear_system_tests.push_back( new SpedicatoFunction17( 500 ) );

    nonlinear_system_tests.push_back( new SSTnonlinearityTerm( 0 ) );
    nonlinear_system_tests.push_back( new SSTnonlinearityTerm( 1 ) );

    nonlinear_system_tests.push_back( new StrictlyConvexFunction1( 2 ) );
    nonlinear_system_tests.push_back( new StrictlyConvexFunction1( 10 ) );
    nonlinear_system_tests.push_back( new StrictlyConvexFunction1( 50 ) );
    nonlinear_system_tests.push_back( new StrictlyConvexFunction1( 500 ) );
    // nonlinear_system_tests.push_back(new StrictlyConvexFunction1(5000));

    nonlinear_system_tests.push_back( new StrictlyConvexFunction2( 2 ) );
    nonlinear_system_tests.push_back( new StrictlyConvexFunction2( 10 ) );
    nonlinear_system_tests.push_back( new StrictlyConvexFunction2( 50 ) );
    nonlinear_system_tests.push_back( new StrictlyConvexFunction2( 500 ) );
    // nonlinear_system_tests.push_back(new StrictlyConvexFunction2(5000));

    nonlinear_system_tests.push_back( new Toint225( 10 ) );
    nonlinear_system_tests.push_back( new Toint225( 100 ) );
    nonlinear_system_tests.push_back( new Toint225( 500 ) );

    nonlinear_system_tests.push_back( new TridimensionalValley() );

    nonlinear_system_tests.push_back( new TrigonometricExponentialSystem1( 10 ) );
    nonlinear_system_tests.push_back( new TrigonometricExponentialSystem1( 50 ) );
    nonlinear_system_tests.push_back( new TrigonometricExponentialSystem1( 500 ) );

    nonlinear_system_tests.push_back( new TrigonometricExponentialSystem2( 9 ) );
    nonlinear_system_tests.push_back( new TrigonometricExponentialSystem2( 27 ) );
    nonlinear_system_tests.push_back( new TrigonometricExponentialSystem2( 81 ) );

    nonlinear_system_tests.push_back( new TrigExp( 10 ) );
    nonlinear_system_tests.push_back( new TrigExp( 100 ) );

    nonlinear_system_tests.push_back( new TrigonometricFunction( 2 ) );
    nonlinear_system_tests.push_back( new TrigonometricFunction( 10 ) );
    nonlinear_system_tests.push_back( new TrigonometricFunction( 50 ) );

    nonlinear_system_tests.push_back( new TroeschFunction( 2 ) );
    nonlinear_system_tests.push_back( new TroeschFunction( 10 ) );
    nonlinear_system_tests.push_back( new TroeschFunction( 50 ) );
    nonlinear_system_tests.push_back( new TroeschFunction( 500 ) );
    // nonlinear_system_tests.push_back(new TroeschFunction(5000));

    nonlinear_system_tests.push_back( new TwoPointBoundaryValueProblem( 10 ) );
    nonlinear_system_tests.push_back( new TwoPointBoundaryValueProblem( 100 ) );
    nonlinear_system_tests.push_back( new TwoPointBoundaryValueProblem( 1000 ) );
    // nonlinear_system_tests.push_back(new TwoPointBoundaryValueProblem(5000));

    nonlinear_system_tests.push_back( new VariablyDimensionedFunction( 5 ) );
    nonlinear_system_tests.push_back( new VariablyDimensionedFunction( 10 ) );
    nonlinear_system_tests.push_back( new VariablyDimensionedFunction( 50 ) );

    nonlinear_system_tests.push_back( new WatsonFunction() );
    nonlinear_system_tests.push_back( new Weibull() );

    nonlinear_system_tests.push_back( new WoodFunction() );

    nonlinear_system_tests.push_back( new XiaoYin1() );
    nonlinear_system_tests.push_back( new XiaoYin2() );
    nonlinear_system_tests.push_back( new XiaoYin3() );

    nonlinear_system_tests.push_back( new YixunShi1() );
    nonlinear_system_tests.push_back( new YixunShi2() );
    nonlinear_system_tests.push_back( new YixunShi3() );
    nonlinear_system_tests.push_back( new YixunShi4() );

    nonlinear_system_tests.push_back( new ZeroJacobianFunction( 10 ) );
    nonlinear_system_tests.push_back( new ZeroJacobianFunction( 50 ) );
    nonlinear_system_tests.push_back( new ZeroJacobianFunction( 101 ) );

    for ( unsigned i{ 0 }; i < unsigned( nonlinear_system_tests_map.size() ); ++i )
    {
      nonlinear_system_tests_map[nonlinear_system_tests[i]->title()] = i;
    }
  }

}  // namespace Utils
