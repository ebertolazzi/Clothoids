%======================================================================%
% FresnelCS:  Compute Fresnel sine and cosine integrals                %
%                                                                      %
% USAGE: [FresnelC,FresnelS] = FresnelCS(y) ;                          %
%                                                                      %
%  Fresnel integral are defined as:                                    %
%                                                                      %
%  C(y) = int_0^y cos( (pi/2)*t^2 ) dt                                 %
%  S(y) = int_0^y sin( (pi/2)*t^2 ) dt                                 %
%                                                                      %
%  The algorithm is described in:                                      %
%    Atlas for computing mathematical functions: an illustrated guide  %
%    for practitioners, with programs in C and Mathematica.            %
%    William J. Thompson New York : Wiley, c1997.                      %
%                                                                      %
%  The code is a sligly modification and translation to C language     %
%  of original code developed by                                       %
%  Venkata Sivakanth Telasula (sivakanth.telasula@gmail.com)           %
%                                                                      %
%  On input:                                                           %
%                                                                      %
%   y = argument of the function C(y) and S(y), it may be a vector     %
%                                                                      %
%  On output:                                                          %
%                                                                      %
%  FresnelC = The value(s) of C(y)                                     %
%  FresnelS = The value(s) of S(y)                                     %
%                                                                      %
%======================================================================%
%                                                                      %
%  Autor: Enrico Bertolazzi                                            %
%         Department of Industrial Engineering                         %
%         University of Trento                                         %
%         enrico.bertolazzi@unitn.it                                   %
%                                                                      %
%======================================================================%
function [FresnelC,FresnelS] = FresnelCS(y)
  error('this function is mex only. Run CompileLib.m script to build')
end