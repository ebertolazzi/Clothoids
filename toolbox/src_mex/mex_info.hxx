/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2016
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#define MEX_INFO_MESSAGE(CMD) \
"\n" \
"  " CMD "( 'delete', OBJ );\n" \
"\n" \
"  " CMD "( 'copy', OBJ, OBJ1 );\n" \
"\n" \
"    [ xmin, ymin, xmax, ymax ] = " CMD "('bbox', varargin );\n" \
"  " CMD "( 'changeOrigin', OBJ, x0, y0 );\n" \
"  " CMD "( 'translate', OBJ, tx, ty );\n" \
"  " CMD "( 'trim', OBJ, smin, smax );\n" \
"  " CMD "( 'rotate', OBJ, angle, cx, cy );\n" \
"  " CMD "( 'scale', OBJ, scale );\n" \
"  " CMD "( 'reverse', OBJ );\n" \
"  " CMD "( 'changeOrigin', OBJ, newX0, newY0 );\n" \
"\n" \
"  [x,y,theta,kappa] = " CMD "( 'evaluate', OBJ, ss );\n" \
"  [X,Y] = " CMD "( 'eval', OBJ, s [,t] );\n" \
"  [X,Y] = " CMD "( 'eval_D', OBJ, s [,t] );\n" \
"  [X,Y] = " CMD "( 'eval_DD', OBJ, s [,t] );\n" \
"  [X,Y] = " CMD "( 'eval_DDD', OBJ, s [,t] );\n" \
"  res = " CMD "( 'theta', OBJ, s );\n" \
"  res = " CMD "( 'theta_D', OBJ, s );\n" \
"  res = " CMD "( 'theta_DD', OBJ, s );\n" \
"  res = " CMD "( 'theta_DDD', OBJ, s );\n" \
"  res = " CMD "( 'kappa', OBJ, s );\n" \
"  res = " CMD "( 'kappa_D', OBJ, s );\n" \
"  res = " CMD "( 'kappa_DD', OBJ, s );\n" \
"\n" \
"  res   = " CMD "( 'xBegin', OBJ );\n" \
"  res   = " CMD "( 'yBegin', OBJ );\n" \
"  [x,y] = " CMD "( 'xyBegin', OBJ );\n" \
"  res   = " CMD "( 'thetaBegin', OBJ );\n" \
"  res   = " CMD "( 'kappaBegin', OBJ );\n" \
"  res   = " CMD "( 'xEnd', OBJ );\n" \
"  res   = " CMD "( 'yEnd', OBJ );\n" \
"  [x,y] = " CMD "( 'xyEnd', OBJ );\n" \
"  res   = " CMD "( 'thetaEnd', OBJ );\n" \
"  res   = " CMD "( 'kappaEnd', OBJ );\n" \
"  res   = " CMD "( 'length', OBJ [,offs] );\n" \
"\n" \
"\n" \
"  [x, y, s, t, iflag, dst] = " CMD "( 'closestPoint', OBJ, x, y [,offs] );\n" \
"  dst = " CMD "( 'distance', OBJ, x, y [,offs] );\n" \
"  [s,t] = " CMD "( 'findST', OBJ, x, y );\n" \
"\n" \
"  ok      = " CMD "( 'collision', OBJ1, OBJ2, [,offs1, offs2] );\n" \
"  [s1,s2] = " CMD "( 'intersect', OBJ1, OBJ2, [,offs1, offs2] );\n" \
"\n" \
"  [s,t] = " CMD "( 'findST', OBJ, x, y );\n" \
"  " CMD "( 'info', OBJ );\n" \
"\n"

#define MEX_INFO_MESSAGE_END \
"=====================================================================================\n" \
"\n" \
"Autors: Enrico Bertolazzi^(1), Marco Frego^(2)\n" \
"  (1) Department of Industrial Engineering\n" \
"  (2) Department of Information Engineering and Computer Science\n" \
"  University of Trento\n" \
"  enrico.bertolazzi@unitn.it\n" \
"  m.fregox@gmail.com\n" \
"\n" \
"=====================================================================================\n"
