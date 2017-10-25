%=============================================================================%
%  buildClothoid:  Compute parameters of the G2 Hermite clothoid fitting      %
%                                                                             %
%  USAGE: [ S0,S1,SM,SG,iter ] = buildClothoid3arcG2( x0, y0, th0, k0,        %
%                                                     x1, y1, th1, k1 ) ;     %
%                                                                             %
%  USAGE: [ S0,S1,SM,SG,iter ] = buildClothoid3arcG2( x0, y0, th0, k0,        %
%                                                     x1, y1, th1, k1,        %
%                                                     dmax0, dmax1 ) ;        %
%                                                                             %
%  On input:                                                                  %
%                                                                             %
%       x0, y0  = coodinate of initial point                                  %
%       theta0  = orientation (angle) of the clothoid at initial point        %
%       k0      = initial curvature                                           %
%       x1, y1  = coodinate of final point                                    %
%       theta1  = orientation (angle) of the clothoid at final point          %
%                                                                             %
%       dmax0   = rough desidered maximum angle divergence from guess         %
%       dmax1   = rough desidered maximum angle variation                     %
%                                                                             %
%  On output:                                                                 %
%                                                                             %
%       S0     = initial arc of clothoid                                      %
%       SM     = middle arc of clothoid                                       %
%       S1     = final arc of clothoid                                        %
%       SG     = G1 guess clothoid                                            %
%                                                                             %
%  Optional Output                                                            %
%                                                                             %
%       iter   = >0 number of iteration used for the computation              %
%                -1 computation failed                                        %
%                                                                             %
%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi and Marco Frego                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          marco.frego@unitn.it                                               %
%                                                                             %
%=============================================================================%
function [ S0,S1,SM,SG,iter ] = buildClothoid3arcG2( x0, y0, th0, k0, x1, y1, th1, k1 )
  error('this function is mex only. Run CompileLib.m script to build')
end
