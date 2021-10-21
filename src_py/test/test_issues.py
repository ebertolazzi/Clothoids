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

RESULT_THETA = [2.1407711050824485, 2.1527204693839153, 2.162658198926908]
RESULT_KAPPA = [0.1295518168070438, 0.10943546922229852, 0.0893191216375532]
RESULT_X = [-496.33422815300224, -496.3887050197814, -496.4440964109352]
RESULT_Y = [2015.865256473188, 2015.9491144003082, 2016.0323713005112]


class TestIssue2(unittest.TestCase):

    def test_evaluate(self):
        XS = [-496.280842337990, -497.213911531027, -497.848398582311, -498.193634183612,
              -498.230956951389, -497.391194677679, -495.319781069411]
        YS = [2015.78070002887, 2017.18363315053, 2018.72410899866, 2020.44339040015,
              2022.39649448451, 2027.22423957009, 2033.63372864574]
        CHALLENGE_S = [0.1, 0.2, 0.3]
        RESULT_THETA = [2.1407711050824485, 2.1527204693839153, 2.162658198926908]
        RESULT_KAPPA = [0.1295518168070438, 0.10943546922229852, 0.0893191216375532]
        RESULT_XS = [-496.33422815300224, -496.3887050197814, -496.4440964109352]
        RESULT_YS = [2015.865256473188, 2015.9491144003082, 2016.0323713005112]

        cl = G2lib.buildP5(XS, YS)
        # Single value
        th, kp, xs, ys = cl.evaluate(CHALLENGE_S[0])
        self.assertAlmostEqual(th, RESULT_THETA[0], places=6)
        self.assertAlmostEqual(kp, RESULT_KAPPA[0], places=6)
        self.assertAlmostEqual(xs, RESULT_XS[0], places=6)
        self.assertAlmostEqual(ys, RESULT_YS[0], places=6)
        # Vector value
        th, kp, xs, ys = cl.evaluate(CHALLENGE_S)
        self.assertEqual(len(RESULT_THETA), len(CHALLENGE_S))
        self.assertEqual(len(RESULT_KAPPA), len(CHALLENGE_S))
        self.assertEqual(len(RESULT_XS), len(CHALLENGE_S))
        self.assertEqual(len(RESULT_YS), len(CHALLENGE_S))
        [self.assertAlmostEqual(a, b, places=6) for a, b in zip(th, RESULT_THETA)]
        [self.assertAlmostEqual(a, b,  places=6) for a, b in zip(kp, RESULT_KAPPA)]
        [self.assertAlmostEqual(a, b,  places=6) for a, b in zip(xs, RESULT_XS)]
        [self.assertAlmostEqual(a, b,  places=6) for a, b in zip(ys, RESULT_YS)]
    
    def test_eval(self):
        XS = [-496.280842337990, -497.213911531027, -497.848398582311, -498.193634183612,
              -498.230956951389, -497.391194677679, -495.319781069411]
        YS = [2015.78070002887, 2017.18363315053, 2018.72410899866, 2020.44339040015,
              2022.39649448451, 2027.22423957009, 2033.63372864574]
        CHALLENGE_S = [0.1, 0.2, 0.3]
        RESULT_XS = [-496.33422815300224, -496.3887050197814, -496.4440964109352]
        RESULT_YS = [2015.865256473188, 2015.9491144003082, 2016.0323713005112]

        cl = G2lib.buildP5(XS, YS)
        # Single value
        xs, ys = cl.eval(CHALLENGE_S[0])
        self.assertAlmostEqual(xs, RESULT_XS[0], places=6)
        self.assertAlmostEqual(ys, RESULT_YS[0], places=6)
        # Vector value
        xs, ys = cl.eval(CHALLENGE_S)
        self.assertEqual(len(xs), len(CHALLENGE_S))
        self.assertEqual(len(ys), len(CHALLENGE_S))
        [self.assertAlmostEqual(a, b,  places=6) for a, b in zip(xs, RESULT_XS)]
        [self.assertAlmostEqual(a, b,  places=6) for a, b in zip(ys, RESULT_YS)]

    def test_findST1(self):
        XS = [-496.280842337990, -497.213911531027, -497.848398582311, -498.193634183612,
              -498.230956951389, -497.391194677679, -495.319781069411]
        YS = [2015.78070002887, 2017.18363315053, 2018.72410899866, 2020.44339040015,
              2022.39649448451, 2027.22423957009, 2033.63372864574]
        CHALLENGE_XS = [-497.0, -496.0, -495.0]
        CHALLENGE_YS = [2015.0, 2016.0, 2022.0]
        RESULT_S = [0.0, 0.03616461826466033, 6.7872905283877945]
        RESULT_N = [0.0, -0.35438664932566455, -3.2455947673245884]

        cl = G2lib.buildP5(XS, YS)
        # Single value
        _, ss, ns = cl.findST1(CHALLENGE_XS[0], CHALLENGE_YS[0])
        self.assertAlmostEqual(ss, RESULT_S[0], places=6)
        self.assertAlmostEqual(ns, RESULT_N[0], places=6)
        # Vector value
        _, ss, ns = cl.findST1(CHALLENGE_XS, CHALLENGE_YS)
        self.assertEqual(len(ss), len(CHALLENGE_XS))
        self.assertEqual(len(ns), len(CHALLENGE_XS))
        [self.assertAlmostEqual(a, b,  places=6) for a, b in zip(ss, RESULT_S)]
        [self.assertAlmostEqual(a, b,  places=6) for a, b in zip(ns, RESULT_N)]


class TestIssue3(unittest.TestCase):

    def setUp(self):
        X = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        Y = [2.0, 1.0, 5.0, 3.0, 4.0, 0.0]
        self.cl = G2lib.buildP5(X, Y)
        self.x = 1.5916858361759525
        self.y = 2.769657501900653
        self.s = 3.8487782381907865
        self.t = -0.46880499443796325
        self.dst = 0.46880499443796336
        self.i = 1

    def test_closestPointInSRange_ISO(self):
        _, x, y, s, t, dst, i = self.cl.closestPointInSRange_ISO(2.0, 3.0, 0.0, self.cl.length() / 2.0)
        self.assertAlmostEqual(x, self.x, places=6)
        self.assertAlmostEqual(y, self.y, places=6)
        self.assertAlmostEqual(s, self.s, places=6)
        self.assertAlmostEqual(t, self.t, places=6)
        self.assertAlmostEqual(dst, self.dst, places=6)
        self.assertEqual(i, self.i)

    def test_closestPointInSRange_SAE(self):
        _, x, y, s, t, dst, i = self.cl.closestPointInSRange_SAE(2.0, 3.0, 0.0, self.cl.length() / 2.0)
        self.assertAlmostEqual(x, self.x, places=6)
        self.assertAlmostEqual(y, self.y, places=6)
        self.assertAlmostEqual(s, self.s, places=6)
        self.assertAlmostEqual(t, -self.t, places=6)
        self.assertAlmostEqual(dst, self.dst, places=6)
        self.assertEqual(i, self.i)


class TestIssue3bis(unittest.TestCase):

    def setUp(self):
        line = G2lib.LineSegment(-1.0, 0.0, 0.0, 2.0)
        self.cl = G2lib.ClothoidList()
        self.cl.push_back(line)

    def test_closestPointInSRange_ISO(self):
        r, x, y, s, t, dst, i = self.cl.closestPointInSRange_ISO(-0.5, 1.0, 0.0, self.cl.length())
        self.assertAlmostEqual(x, -0.5, places=6)
        self.assertAlmostEqual(y, 0.0, places=6)
        self.assertAlmostEqual(s, 0.5, places=6)
        self.assertAlmostEqual(t, 1.0, places=6)
        self.assertAlmostEqual(dst, 1.0, places=6)
        self.assertEqual(i, 0)
        self.assertEqual(r, 1)


class TestIssue4(unittest.TestCase):

    def test_eval_ISO(self):
        line = G2lib.LineSegment(0.0, 0.0, 0.0, 0.0)
        line.build_2P([0, 0], [10, 0])
        SS = [line.length() / 5.0 * i for i in range(6)]
        NS = [i for i in range(6)]
        US = [2*i for i in range(6)]
        VS = NS
        xs, ys = line.eval_ISO(SS, NS)
        for x, u in zip(xs, US):
            self.assertAlmostEqual(x, u, places=5)
        for y, v in zip(ys, VS):
            self.assertAlmostEqual(y, v, places=5)
            

class TestIssue5(unittest.TestCase):

    def setUp(self):
        self.x = [355.657984270481, 350.227521567838, 348.286737646442, 346.812488321448, 
             345.786112209316, 345.160955850035, 344.918357859831, 345.030326162931]
        self.y = [-245.287803483196, -245.397808791138, -245.109044851735, -244.586519603618, 
             -243.802731677890, -242.743930327706, -241.396364783868, -239.732533575036]
        self.cloth = G2lib.buildP5(self.x, self.y)

    def test_getSK(self):
        s = [0.0, 5.434222626968327, 7.398631315590526, 8.966846624943887, 10.2653773327908, 
             11.502074712176352, 12.878055135384317, 14.546315256266153]
        k = [0.056762709479033974, -0.07666886061868458, -0.09249440048977871, -0.2254123988506045, -0.3326350963243614,
             -0.26845631817438453, -0.23062166055358554, 0.18079726534322924]
        s_, k_ = self.cloth.getSK()
        for u, v in zip(s_, s):
            self.assertAlmostEqual(u, v, places=5)
        for u, v in zip(k_, k):
            self.assertAlmostEqual(u, v, places=5)

    def test_getSTK(self):
        s = [0.0, 5.434222626968327, 7.398631315590526, 8.966846624943887, 10.2653773327908, 
             11.502074712176352, 12.878055135384317, 14.546315256266153]
        k = [0.056762709479033974, -0.07666886061868458, -0.09249440048977871, -0.2254123988506045, -0.3326350963243614,
             -0.26845631817438453, -0.23062166055358554, 0.18079726534322924]
        t = [-3.1547245175642336, -3.2088117460336982, -3.3749646359922596, -3.6242377908288153, -3.9865586952896903, 
             -4.358242783830572, -4.701603548022518, -4.743163573810764]
        s_, t_, k_ = self.cloth.getSTK()
        for u, v in zip(s_, s):
            self.assertAlmostEqual(u, v, places=5)
        for u, v in zip(t_, t):
            self.assertAlmostEqual(u, v, places=5)
        for u, v in zip(k_, k):
            self.assertAlmostEqual(u, v, places=5)
    
    def test_getXY(self):
        x_, y_ = self.cloth.getXY()
        for u, v in zip(x_, self.x):
            self.assertAlmostEqual(u, v, places=5)
        for u, v in zip(y_, self.y):
            self.assertAlmostEqual(u, v, places=5)

    def test_getDeltaTheta(self):
        dtheta = [0.0, -8.881784197001252e-16, 8.881784197001252e-16, 
                  4.440892098500626e-16, 0.0, 0.0]
        dtheta_ = self.cloth.getDeltaTheta()
        for u, v in zip(dtheta_, dtheta):
           self.assertAlmostEqual(u, v, places=5)
    
    def test_getDeltaKappa(self):
        dkappa = [-2.0816681711721685e-15, -5.093148125467906e-15, -1.887379141862766e-15, 
                  4.107825191113079e-15, 7.771561172376096e-15, 2.7755575615628914e-15]
        dkappa_ = self.cloth.getDeltaKappa()
        for u, v in zip(dkappa_, dkappa):
           self.assertAlmostEqual(u, v, places=5)


# class TestIssue6(unittest.TestCase):
# 
#     def test_findST_ISO(self):
#         import json
#         from os import path
# 
#         with open(path.normpath(path.join(__file__, "..", "data_test_issue_6.json"))) as dp:
#             data = json.load(dp)
# 
#         cloth = G2lib.buildP1(data["x_cl"], data["y_cl"], data["m0_cl"], data["m1_cl"])
#         ss, _ = cloth.findST_ISO(data["x_query"], data["y_query"])
# 
#         query_size = len(ss)
#         self.assertEqual(query_size, len(data["x_query"]))
#         for i in range(1, query_size):
#             self.assertTrue(ss[i] - ss[i-1] >= 0)



