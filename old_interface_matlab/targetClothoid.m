%======================================================================%
%  bbClothoid:  Compute total variation, curvature and jerk            %
%                                                                      %
%  USAGE:                                                              %
%    [tv,curv2,jerk2] = targetClothoid( S ) ;                          %
%    [tv,curv2,jerk2] = targetClothoid( x0, y0, theta0, k, dk, L ) ;   %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%    x0, y0 = coodinate of initial point                               %
%    theta0 = orientation (angle) of the clothoid at initial point     %
%    k      = curvature at initial point                               %
%    dk     = derivative of curvature respect to arclength             %
%    L      = the lenght of the clothoid curve or a vector of length   %
%             where to compute the clothoid values                     %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function [tv,curv2,jerk2] = targetClothoid( S )
  error('this function is mex only. Run CompileLib.m script to build')
end