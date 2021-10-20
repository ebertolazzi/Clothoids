# PYTHON Wrapper for Clothoids
# 
# License MIT - See LICENSE file
# 
# 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
#      Enrico Bertolazzi, Marco Frego

import sys
import os
if sys.platform != 'win32':
    sys.path.insert(0, os.path.normpath(os.path.join(__file__, "../../build")))
else:
    sys.path.insert(0, os.path.normpath(os.path.join(__file__, "../../build/Release")))

import unittest
import G2lib

XS = [
    -496.280842337990, -497.213911531027, -497.848398582311, -498.193634183612,
    -498.230956951389, -497.391194677679, -495.319781069411 
]

YS = [
    2015.78070002887, 2017.18363315053, 2018.72410899866, 2020.44339040015,
    2022.39649448451, 2027.22423957009, 2033.63372864574 
]

THETA_0 = 2.24459953434543
THETA_1 = 1.21583530295506

if sys.platform != 'win32':
    class TestInterpolator(unittest.TestCase):

        @staticmethod
        def get_theta(clothlist):
            return [cc.thetaBegin() for cc in clothlist.as_list()]


        def test_P1(self):
            THETA_P1 = [2.2446, 2.0618, 1.8645, 1.6764, 1.5165, 1.3143, 1.2158]
            clothlist = G2lib.buildP1(XS, YS, THETA_0, THETA_1)
            theta_opt = self.get_theta(clothlist)
            for t_opt, t_test in zip(theta_opt, THETA_P1):
                self.assertAlmostEqual(t_opt, t_test, places=4,
                msg="Difference grater than 1e-4 in theta for P1")
        
        def test_P2(self):
            THETA_P2 = [2.0088, 2.1245, 1.8488, 1.6748, 1.5409, 1.1464, 2.0088]
            clothlist = G2lib.buildP2(XS, YS)
            theta_opt = self.get_theta(clothlist)
            for t_opt, t_test in zip(theta_opt, THETA_P2):
                self.assertAlmostEqual(t_opt, t_test, places=4,
                msg="Difference grater than 1e-4 in theta for P2")

        def test_P4(self):
            THETA_P4 = [2.2568, 2.0586, 1.8654, 1.6763, 1.5160, 1.3181, 1.1983]
            clothlist = G2lib.buildP4(XS, YS)
            theta_opt = self.get_theta(clothlist)
            for t_opt, t_test in zip(theta_opt, THETA_P4):
                self.assertAlmostEqual(t_opt, t_test, places=4,
                msg="Difference grater than 1e-4 in theta for P4")
        
        def test_P5(self):
            THETA_P5 = [2.1268, 2.0934, 1.8559, 1.6784, 1.5177, 1.3009, 1.2773]
            clothlist = G2lib.buildP5(XS, YS)
            theta_opt = self.get_theta(clothlist)
            for t_opt, t_test in zip(theta_opt, THETA_P5):
                self.assertAlmostEqual(t_opt, t_test, places=4,
                msg="Difference grater than 1e-4 in theta for P5")

        def test_P6(self):
            THETA_P6 = [2.1581, 2.0850, 1.8582, 1.6779, 1.5172, 1.3052, 1.2579]
            clothlist = G2lib.buildP6(XS, YS)
            theta_opt = self.get_theta(clothlist)
            for t_opt, t_test in zip(theta_opt, THETA_P6):
                self.assertAlmostEqual(t_opt, t_test, places=4,
                msg="Difference grater than 1e-4 in theta for P6")

        def test_P7(self):
            THETA_P7 = [2.1996, 2.0739, 1.8612, 1.6773, 1.5167, 1.3108, 1.2319]
            clothlist = G2lib.buildP7(XS, YS)
            theta_opt = self.get_theta(clothlist)
            for t_opt, t_test in zip(theta_opt, THETA_P7):
                self.assertAlmostEqual(t_opt, t_test, places=4,
                msg="Difference grater than 1e-4 in theta for P7")

        def test_P8(self):
            THETA_P8 = [2.2570, 2.0585, 1.8654, 1.6763, 1.5162, 1.3163, 1.2064]
            clothlist = G2lib.buildP8(XS, YS)
            theta_opt = self.get_theta(clothlist)
            for t_opt, t_test in zip(theta_opt, THETA_P8):
                self.assertAlmostEqual(t_opt, t_test, places=4,
                msg="Difference grater than 1e-4 in theta for P8")

        def test_P9(self):
            THETA_P9 = [2.2563, 2.0587, 1.8653, 1.6763, 1.5162, 1.3163, 1.2065]
            clothlist = G2lib.buildP9(XS, YS)
            theta_opt = self.get_theta(clothlist)
            for t_opt, t_test in zip(theta_opt, THETA_P9):
                self.assertAlmostEqual(t_opt, t_test, places=4,
                msg="Difference grater than 1e-4 in theta for P9")


    