%======================================================================%
%  pointsOnClothoid:  Compute points on a clothoid curve.              %
%                     Used for plotting purpose.                       %
%                                                                      %
%  USAGE:                                                              %
%    XY    = pointsOnClothoid( x0, y0, theta0, k0, dk, ss ) ;          %
%    [X,Y] = pointsOnClothoid( x0, y0, theta0, k0, dk, ss ) ;          %
%    XY    = pointsOnClothoid( x0, y0, theta0, k0, dk, ss, offs) ;     %
%    [X,Y] = pointsOnClothoid( x0, y0, theta0, k0, dk, ss, offs) ;     %
%    XY    = pointsOnClothoid( clot, ss ) ;                            %
%    [X,Y] = pointsOnClothoid( clot, ss ) ;                            %
%    XY    = pointsOnClothoid( clot, ss, offs ) ;                      %
%    [X,Y] = pointsOnClothoid( clot, ss, offs ) ;                      %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%    x0, y0  = coodinate of initial point                              %
%    theta0  = orientation (angle) of the clothoid at initial point    %
%    k0      = curvature at initial point                              %
%    dk      = derivative of curvature respect to arclength            %
%    ss      = a vector with the curvilinear coordinates where         %
%              to compute the clothoid values                          %
%    offs    = curve offset                                            %
%                                                                      %
%    clot    = struct with field `x0`, `y0`, `theta0`, `k0`, `dk`      %
%              and a final field                                       %
%                `L` the lenght of the clothoid curve                  %
%              or in alternative two field                             %
%                 `smin` = initial curvilinear coordinate of the curve %
%                 `smax` = final curvilinear coordinate of the curve   %
%                                                                      %
%  On output: (1 argument)                                             %
%                                                                      %
%    XY = matrix 2 x NPTS whose column are the points of the clothoid  %
%                                                                      %
%  On output: (2 argument)                                             %
%                                                                      %
%    X  = matrix 1 x NPTS X coordinate of points of the clothoid       %
%    Y  = matrix 1 x NPTS Y coordinate of points of the clothoid       %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function XY    = pointsOnClothoid( x0, y0, theta0, k0, dk, ss )
  error('this function is mex only. Run CompileLib.m script to build')
end