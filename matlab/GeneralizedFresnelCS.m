%===================================================================%
%  Compute Fresnel sine and cosine integrals momenta                %
%                                                                   %
%  USAGE:                                                           %
%    [X,Y] = GeneralizedFresnelCS( nk, a, b, c ) ;                  %
%                                                                   %
%  Integrals are defined as:                                        %
%                                                                   %
%    X_k(a,b,c) = int_0^1 t^k * cos( (a/2)*t^2 + b*t + c ) dt       %
%    Y_k(a,b,c) = int_0^1 t^k * sin( (a/2)*t^2 + b*t + c ) dt       %
%                                                                   %
%  On input:                                                        %
%                                                                   %
%    nk      = number of momentae to be computed (1 <= nk <= 3)     %
%    a, b, c = the parameters of the integrals or vectors of the    %
%              same dimensions length(a)=length(b)=length(c)        %
%                                                                   %
%  On output:                                                       %
%                                                                   %
%    X = vector with Fresnel cosine momenta [X_0,X_1,...,X_{nk-1}]  %
%    Y = vector with Fresnel sine momenta   [Y_0,Y_1,...,Y_{nk-1}]  %
%                                                                   %
%    the dimension of the vectors X are nk x length(a)              %
%                                                                   %
%===================================================================%
%                                                                   %
%  Autor: Enrico Bertolazzi                                         %
%         Department of Industrial Engineering                      %
%         University of Trento                                      %
%         enrico.bertolazzi@unitn.it                                %
%                                                                   %
%===================================================================%
function [X,Y] = GeneralizedFresnelCS( nk, a, b, c )
  error('this function is mex only. Run CompileLib.m script to build')
end

