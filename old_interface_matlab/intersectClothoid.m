%======================================================================%
%  intersectClothoid:  Compute intersections betweed clothoids         %
%                                                                      %
%  USAGE:                                                              %
%    [s1, s2] = intersectClothoid( clot1, clot2 ) ;                    %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%    clot1 = structure with the definition of the first clothoid       %
%    clot2 = structure with the definition of the second clothoid      %
%                                                                      %
%    The structure must contain the field                              %
%      x0     = initial x coordinate                                   %
%      y0     = initial y coordinate                                   %
%      theta0 = orientation (angle) of the clothoid at initial point   %
%      k0     = curvature at initial point                             %
%      dk     = derivative of curvature respect to arclength           %
%                                                                      %
%    A final field                                                     %
%      L      = the lenght of the clothoid curve                       %
%                                                                      %
%    or in alternative                                                 %
%      smin   = initial curvilinear coordinate of the curve            %
%      smax   = final curvilinear coordinate of the curve              %
%                                                                      %
%    NOTE: (x0,y0,theta0,kappa0) are at curviliear coodinates 0        %
%          smax  > smin  and smin may be negative                      %
%          specify L is equivalent to pass smin = 0 and smax = L       %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%  s1 = curvilinear coordinates of intersections on clot1              %
%  s2 = curvilinear coordinates of intersections on clot2              %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function [s1, s2] = intersectClothoid( clot1, clot2 ) ;
  error('this function is mex only. Run CompileLib.m script to build')
end