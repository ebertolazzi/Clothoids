%======================================================================%
%  bbClothoid:  Compute a series of bounding triangles for             %
%               a clothoid curve.                                      %
%                                                                      %
%  USAGE:                                                              %
%    TT = bbClothoid( clot, angle, size ) ;                            %
%    TT = bbClothoid( clot, angle, size, offs ) ;                      %
%    TT = bbClothoid( x0, y0, theta0, k0, dk, L, angle, size ) ;       %
%    TT = bbClothoid( x0, y0, theta0, k0, dk, L, angle, size, offs ) ; %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%  x0, y0 = coodinate of initial point                                 %
%  theta0 = orientation (angle) of the clothoid at initial point       %
%  k0     = curvature at initial point                                 %
%  dk     = derivative of curvature respect to arclength               %
%  L      = the lenght of the clothoid curve or a vector of length     %
%           where to compute the clothoid values                       %
%  angle  = maximum variation of angle in the bounding box triangle    %
%  size   = maximum height of the bounding box triangle                %
%                                                                      %
%  clot    = struct with field `x0`, `y0`, `theta0`, `k0`, `dk`        %
%            and a final field                                         %
%              `L` the lenght of the clothoid curve                    %
%            or in alternative two field                               %
%               `smin` = initial curvilinear coordinate of the curve   %
%               `smax` = final curvilinear coordinate of the curve     %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%  TT = matrix 6 x N whose column are the coordinates of bounding      %
%       triangles.                                                     %
%       TT(:,i) = [ x0, y0, x1, y1, x2, y2 ].'                         %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function TT = bbClothoid( clot, angle, size )
  error('this function is mex only. Run CompileLib.m script to build')
end