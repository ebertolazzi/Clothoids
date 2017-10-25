%=============================================================================%
%  buildClothoid:  Compute parameters of the G1 Hermite clothoid fitting      %
%                                                                             %
%  USAGE:                                                                     %
%    [S,iter] = buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;             %
%    [S,iter] = buildClothoid( x0, y0, theta0, x1, y1, theta1, derivative ) ; %
%                                                                             %
%  On input:                                                                  %
%                                                                             %
%    x0, y0     = coodinate of initial point                                  %
%    theta0     = orientation (angle) of the clothoid at initial point        %
%    x1, y1     = coodinate of final point                                    %
%    theta1     = orientation (angle) of the clothoid at final point          %
%    derivative = if present and true compute additional                      %
%                 partial derivative of the solution                          %
%                                                                             %
%  On output:                                                                 %
%                                                                             %
%    S.x     = x-coodinate of initial point                                   %
%    S.y     = y-coodinate of initial point                                   %
%    S.theta = orientation (angle) of the clothoid at initial point           %
%    S.k     = curvature at initial point                                     %
%    S.dk    = derivative of curvature respect to arclength,                  %
%              notice that curvature at final point is k+dk*L                 %
%    S.L     = the lenght of the clothoid curve from initial to final point   %
%    S.iter  = Newton Iterations used to solve the interpolation problem      %
%                                                                             %
%  optional output                                                            %
%                                                                             %
%    S.k_1   = partial derivative of the solution respect to theta0           %
%    S.dk_1  = partial derivative of the solution respect to theta0           %
%    S.L_1   = partial derivative of the solution respect to theta0           %
%                                                                             %
%    S.k_2   = partial derivative of the solution respect to theta1           %
%    S.dk_2  = partial derivative of the solution respect to theta1           %
%    S.L_2   = partial derivative of the solution respect to theta1           %
%                                                                             %
%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi and Marco Frego                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
function [S,iter] = buildClothoid( x0, y0, theta0, x1, y1, theta1 )
  error('this function is mex only. Run CompileLib.m script to build')
end