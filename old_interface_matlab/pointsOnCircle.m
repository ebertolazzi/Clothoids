%======================================================================%
%  pointsOnCircle:  Compute points on a circle curve.                  %
%                   Used for plotting purpose.                         %
%                                                                      %
%  USAGE: XY    = pointsOnCicle( x0, y0, theta0, k, ss ) ;             %
%  USAGE: [X,Y] = pointsOnCicle( x0, y0, theta0, k, ss ) ;             %
%  USAGE: XY    = pointsOnCicle( S, ss ) ;                             %
%  USAGE: [X,Y] = pointsOnCicle( S, ss ) ;                             %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%    S       = struct with field `x`, `y`, `theta`, `k`, `L`           %
%    x0, y0  = coodinate of initial point                              %
%    theta0  = orientation (angle) of the Circle at initial point      %
%    k       = curvature                                               %
%    ss      = the lenght of the circle curve or a vector of length    %
%              where to compute the circle values                      %
%                                                                      %
%  On output: (1 argument)                                             %
%                                                                      %
%    XY = matrix 2 x NPTS whose column are the points of the circle    %
%                                                                      %
%  On output: (2 argument)                                             %
%                                                                      %
%    X  = matrix 1 x NPTS X coordinate of points of the circle         %
%    Y  = matrix 1 x NPTS Y coordinate of points of the circle         %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function XY = pointsOnCicle( x0, y0, theta0, k, ss )
  error('this function is mex only. Run CompileLib.m script to build')
end