%======================================================================%
%  biarc:  Compute biarc fitting.                                      %
%                                                                      %
%  USAGE:                                                              %
%   [arc1,arc2,ok] = biarc( x0, y0, theta0, x1, y1, theta1 ) ;         %
%   [arc1,arc2,ok] = biarc( x0, y0, theta0, x1, y1, theta1, thstar ) ; %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%  x0, y0 = coordinate of initial point                                %
%  theta0 = orientation (angle) at the initial point                   %
%  x1, y1 = coordinate of final point                                  %
%  theta1 = orientation (angle) at the final point                     %
%  thstar = orientation (angle) at the intermediate (optional)         %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%  arc1, arc2 = rational B-spline of the two arc                       %
%  ok         = false if computation fails                             %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function [arc1,arc2,ok] = biarc( x0, y0, theta0, x1, y1, theta1 )
  error('this function is mex only. Run CompileLib.m script to build')
end