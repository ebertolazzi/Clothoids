%=============================================================================%
%  intXY:  Compute Fresnel sine and cosine integrals momenta                  %
%                                                                             %
%  USAGE: [X,Y] = GeneralizedFresnelCS( nk, a, b, c ) ;                       %
%                                                                             %
%  Integrals are defined as:                                                  %
%                                                                             %
%  X_k(a,b,c) = \int_0^1 t^k * cos( (a/2)*t^2 + b*t + c ) dt                  %
%  Y_k(a,b,c) = \int_0^1 t^k * sin( (a/2)*t^2 + b*t + c ) dt                  %
%                                                                             %
%  On input:                                                                  %
%                                                                             %
%       nk      = number of momentae to be computed                           %
%       a, b, c = the parameters of the integrals                             %
%                                                                             %
%  On output:                                                                 %
%                                                                             %
%       X = vector with Fresnel cosine momenta [X_0,X_1,...,X_{nk-1}]         %
%       Y = vector with Fresnel sine momenta   [Y_0,Y_1,...,Y_{nk-1}]         %
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
function [X,Y] = GeneralizedFresnelCS_new( nk, a, b, c )
  
  k = 0;
  if (abs(b)^6/720) < eps
    k = 1 ; % b - small
  end
  if (abs(a)^6/46080) < eps
    k = k + 2 ; % a - small
  end
  
  switch ( k )
  case 0
    [X,Y] = GeneralizedFresnel_ab_large( nk, a, b ) ;
  case 1
    [X,Y] = GeneralizedFresnel_b_small( nk, a, b ) ;
  case 2
    [X,Y] = GeneralizedFresnel_a_small( nk, a, b ) ;
  case 3
    [X,Y] = GeneralizedFresnel_ab_small( nk, a, b ) ;
    %[X,Y] = GeneralizedFresnel_a_small( nk, a, b ) ;
    %[X,Y] = GeneralizedFresnel_b_small( nk, a, b ) ;
  end

  cc = cos(c) ;
  ss = sin(c) ;

  for k=1:nk
    xx = X(k) ;
    yy = Y(k) ;
    X(k) = xx*cc-yy*ss ;
    Y(k) = xx*ss+yy*cc ;
  end
end
%
%
%
function [S,C] = sincos( n, x )
  S    = zeros(n,1) ;
  C    = zeros(n,1) ;
  x2   = x*x ;
  S(1) = sin(x) ;
  C(1) = cos(x) ;
  if n < 2 ; return ; end
  if x2*x2*x2 < eps*5040
    % approx by taylor series
    S(2) = 1-x2*(1-x2/20)/6 ;
    C(2) = 1-x2*(1-x2/30)/12 ;
    if n < 3 ; return ; end
    S(3) = 1-x2*(1-x2/42)/20 ;
    C(3) = 1-x2*(1-x2/56)/30 ;
    if n < 4 ; return ; end
    S(4) = 1-x2*(1-x2/72)/42 ;
    C(4) = 1-x2*(1-x2/90)/56 ;
    if n < 5 ; return ; end
    S(5) = 1-x2*(1-x2/110)/72 ;
    C(5) = 1-x2*(1-x2/132)/90 ;
    if n < 6 ; return ; end
    S(6) = 1-x2*(1-x2/156)/110 ;
    C(6) = 1-x2*(1-x2/182)/132 ;
    if n < 7 ; return ; end
    S(7) = 1-x2*(1-x2/210)/156 ;
    C(7) = 1-x2*(1-x2/240)/182 ;
    if n < 8 ; return ; end
    S(8) = 1-x2*(1-x2/272)/210 ;
    C(8) = 1-x2*(1-x2/306)/240 ;
    if n < 9 ; return ; end
    S(9) = 1-x2*(1-x2/342)/272 ;
    C(9) = 1-x2*(1-x2/380)/306 ;
  else
    if n < 2 ; return ; end
    S(2) = sin(x)/x ;
    C(2) = 2*(1-C(1))/x2 ;
    if n < 3 ; return ; end
    S(3) = 6*(1-S(2))/x2 ;
    C(3) = 12*(1-C(2))/x2 ;
    if n < 4 ; return ; end
    S(4) = 20*(1-S(3))/x2 ;
    C(4) = 30*(1-C(3))/x2 ;
    if n < 5 ; return ; end
    S(5) = 42*(1-S(4))/x2 ;
    C(5) = 56*(1-C(4))/x2 ;
    if n < 6 ; return ; end
    S(6) = 72*(1-S(5))/x2 ;
    C(6) = 90*(1-C(5))/x2 ;
    if n < 7 ; return ; end
    S(7) = 110*(1-S(6))/x2 ;
    C(7) = 132*(1-C(6))/x2 ;
    if n < 8 ; return ; end
    S(8) = 156*(1-S(7))/x2 ;
    C(8) = 182*(1-C(7))/x2 ;
    if n < 9 ; return ; end
    S(9) = 210*(1-S(8))/x2 ;
    C(9) = 240*(1-C(8))/x2 ;
  end
end
%
%
%
function [X,Y] = int_sincos_tau_k( nk, b )
  sz = [ 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9 ] ;

  [S,C] = sincos( sz(nk), b ) ;

  X = zeros(nk,1) ;
  Y = zeros(nk,1) ;

  X(1) = S(2) ;
  Y(1) = C(2)*b/2 ; % 2
  
  if nk < 2 ; return ; end ;

  X(2) = S(2)-C(2)/2 ;
  Y(2) = (C(2)/2-S(3)/6)*b ; % 3

  if nk < 3 ; return ; end ;
  
  X(3) = S(2)-C(2)+S(3)/3 ;
  Y(3) = (C(2)/2-S(3)/3+C(3)/12)*b ; % 3

  if nk < 4 ; return ; end ;

  X(4) = S(2)-1.5*C(2)+S(3)-0.25*C(3) ;
  Y(4) = (C(2)/2-S(3)/2+C(3)/4-S(4)/20)*b ; % 4

  if nk < 5 ; return ; end ;

  X(5) = S(2)-C(3)+2*(S(3)-C(2))+S(4)/5 ;
  Y(5) = (C(2)/2-(2/3)*S(3)+0.5*C(3)-S(4)/5+C(4)/30)*b ; % 4

  if nk < 6 ; return ; end ;

  X(6) = S(2)-2.5*(C(2)+C(3))+(10/3)*S(3)-C(4)/6+S(4) ;
  Y(6) = (0.5*(C(2)-S(4))+(5/6)*(C(3)-S(3))+C(4)/6-S(5)/42)*b ; % 5

  if nk < 7 ; return ; end ;

  X(7) = S(2)-C(4)+3*(S(4)-C(2))+5*(S(3)-C(3))+S(5)/7 ;
  Y(7) = (0.5*(C(2)+C(4))-S(3)+(5/4)*C(3)-S(4)-S(5)/7+C(5)/56)*b ; % 5

  if nk < 8 ; return ; end ;

  X(8) = S(2)-3.5*(C(2)+C(4))+7*(S(3)+S(4))-8.75*C(3)-C(5)/8+S(5) ;
  Y(8) = (0.5*(C(2)-S(5))+(7/6)*(C(4)-S(3))+(7/4)*(C(3)-S(4))+C(5)/8-S(6)/72)*b ; % 6

  if nk < 9 ; return ; end ;

  X(9) = S(2)+4*(S(5)-C(2))+(28/3)*(S(3)-C(4))+14*(S(4)-C(3))-C(5)+S(6)/9 ;
  Y(9) = (0.5*(C(2)+C(5))-(4/3)*(S(3)+S(5))+(7/3)*(C(3)+C(4))-2.8*S(4)-S(6)/9+C(6)/90)*b ;
  
  if nk < 10 ; return ; end ;

  X(10) = S(2)-4.5*(C(2)+C(5))-21*(C(3)+C(4))+12*(S(3)+S(5))+25.2*S(4)-C(6)/10+S(6) ;
  Y(10) = (0.5*(C(2)-S(6))+1.5*(C(5)-S(3))+3*C(3)+4.2*(C(4)-S(4))-3*S(5)+C(6)/10-S(7)/110)*b ; % 7

  if nk < 11 ; return ; end ;

  X(11) = S(2)+5*(S(6)-C(2))+30*(S(5)-C(3))+15*(S(3)-C(5))+42*(S(4)-C(4))-C(6)+S(7)/11 ;
  Y(11) = (0.5*(C(2)+C(6))-(5/3)*(S(3)+S(6))+(15/4)*(C(3)+C(5))-6*(S(4)+S(5))+7*C(4)-S(7)/11+C(7)/132)*b ;

  if nk < 12 ; return ; end ;

  X(12) = S(2)-5.5*(C(2)+C(6))+(55/3)*(S(3)+S(6))-(165/4)*(C(3)+C(5))-77*C(4)+66*(S(4)+S(5))-C(7)/12+S(7) ;
  Y(12) = (0.5*(C(2)-S(7))+(11/6)*(C(6)-S(3))+(55/12)*C(3)+(33/4)*(C(5)-S(4))+11*C(4)-11*S(5)-(55/12)*S(6)+C(7)/12-S(8)/156)*b ; % 8

  if nk < 13 ; return ; end ;

  X(13) = S(2)-6*C(2)-55*C(3)+22*S(3)-132*C(4)+99*S(4)-99*C(5)+132*S(5)-22*C(6)+55*S(6)-C(7)+6*S(7)+S(8)/13 ;
  Y(13) = (0.5*(C(2)+C(7))-2*S(3)+5.5*(C(3)+C(6))-11*S(4)+(33/2)*(C(4)+C(5))-(132/7)*S(5)-11*S(6)-2*S(7)-S(8)/13+C(8)/182)*b ;

  if nk < 14 ; return ; end ;

  X(14) = S(2)-(13/2)*(C(2)+C(7))+26*S(3)-(143/2)*(C(3)+C(6))-(429/2)*(C(4)+C(5))+143*(S(4)+S(6))+(1716/7)*S(5)+26*S(7)-C(8)/14+S(8) ;
  Y(14) = (0.5*(C(2)-S(8))+(13/6)*(C(7)-S(3))+(13/2)*(C(3)-S(7))+(143/10)*(C(6)-S(4))+(143/6)*(C(4)-S(6))+(429/14)*(C(5)-S(5))+C(8)/14-S(9)/210)*b ;  % 9

  if nk < 15 ; return ; end ;

  X(15) = S(2)-7*C(2)+(91/3)*(S(3)-C(7))+91*(S(7)-C(3))+(1001/3)*(S(6)-C(4))+(1001/5)*(S(4)-C(6))+429*(S(5)-C(5))-C(8)+7*S(8)+S(9)/15 ;
  Y(15) = (0.5*(C(2)+C(8))-(7/3)*(S(3)+S(8))+(91/12)*(C(3)+C(7))-(91/5)*(S(4)+S(7))+(1001/30)*(C(4)+C(6))-(143/3)*(S(5)+S(6))+(429/8)*C(5)-S(9)/15+(1/240)*C(9))*b ;

  if nk > 15
    error( 'nk > 15 ') ;
  end ;

end
%
% Xk = int_0^1 cos( a/2 * tau^2 ) * tau^k dtau
% Yk = int_0^1 sin( a/2 * tau^2 ) * tau^k dtau
%
function [X,Y] = int_sincos2_tau_k( nk, aa, b )

  X = zeros(nk,1) ;
  Y = zeros(nk,1) ;

  a  = abs(aa) ;
  as = sqrt(a/pi) ;
  
  SIGNA = 1 ;
  if a < 0
    SIGNA = -1 ;
  end

  [C,S] = FresnelCS(as);

  c0 = C / as ;
  s0 = S / as ;
  
  X(1) = c0 ;
  Y(1) = SIGNA*s0 ;
  
  if nk < 2 ; return ; end ;

  c1 = sin(a/2)/a;
  s1 = (1-cos(a/2))/a ;
  
  X(2) = c1 ;
  Y(2) = SIGNA*s1 ;
  
  if nk < 3 ; return ; end ;

  c2 = (a * c1 - s0) / a;
  s2 = (a * s1 + c0 - 1)/a;
  
  X(3) = c2 ;
  Y(3) = SIGNA*s2 ;
  
  if nk < 4 ; return ; end ;

  c3 = (a * c1 - 2 * s1) / a;
  s3 = (a * s1 + 2 * c1 - 1) / a;
  
  X(4) = c3 ;
  Y(4) = SIGNA*s3 ;
  
  if nk < 5 ; return ; end ;

  c4 = (a * c1 + 3*(1-a * s1 - c0) / a) / a;
  s4 = (a * s1 - 1 + 3*(a * c1 - s0) / a) / a;
  
  X(5) = c4 ;
  Y(5) = SIGNA*s4 ;

  if nk < 6 ; return ; end ;

  c5 = (a * c1 + 4*(1 - a * s1 - 2 * c1) / a) / a;
  s5 = (a * s1 - 1 + 4*( a * c1 - 2 * s1) / a) / a;
  
  X(6) = c5 ;
  Y(6) = SIGNA*s5 ;

  if nk < 7 ; return ; end ;

  c6 = (a * c1 + (-5 * a * s1 + 5 + 15*(s0 - a * c1) / a) / a) / a;
  s6 = (a * s1 - 1 + 5*(a * c1 + 5*(1 - a * s1 - c0) / a) / a) / a;

  X(7) = c6 ;
  Y(7) = SIGNA*s6 ;

  if nk < 8 ; return ; end ;

  c7 = (a * c1 + 6*(1-a * s1 + 4*(2* s1 - a * c1) / a) / a) / a;
  s7 = (a * s1 - 1 + 6*(a * c1 + 4*(1-a * s1 - 2* c1) / a) / a) / a;

  X(8) = c7 ;
  Y(8) = SIGNA*s7 ;

  if nk < 9 ; return ; end ;

  c8 = (a * c1 + 7*(1-a * s1 + 5*(3*(a * s1 + c0 -1) / a - a * c1) / a) / a) / a;
  s8 = (a * s1 - 1 + (7 * a * c1 + 35*(1-a * s1 + 3*(s0 - a * c1) / a) / a) / a) / a;

  X(9) = c8 ;
  Y(9) = SIGNA*s8 ;

  if nk < 10 ; return ; end ;

  c9 = (a * c1 + 8*(1-a * s1 + 6*(4*(a * s1 + 2 * c1 - 1) / a - a * c1) / a) / a) / a;
  s9 = (a * s1 - 1 + 8*(a * c1 + 6*(1-a * s1 + 4*(2* s1-a * c1) / a) / a) / a) / a;

  X(10) = c9 ;
  Y(10) = SIGNA*s9 ;

  if nk > 10
    error( 'nk > 10 ') ;
  end ;

end
%
%
%
function [X,Y] = GeneralizedFresnel_a_small( nk, a, b )

  X = zeros(nk,1) ;
  Y = zeros(nk,1) ;

  [Cb,Sb] = int_sincos_tau_k( 10+nk, b ) ;

  % cos(a/2*tau^2) ~ 1-(1/8)*a^2*tau^4+(1/384)*a^4*tau^8 ;
  % sin(a/2*tau^2) ~  (1/2)*a*tau^2-(1/48)*a^3*tau^6+(1/3840)*a^5*tau^10
  a2 = a*a ;
  for kk=0:nk-1
    Cc = Cb(1+kk) - (a2/8)*( Cb(5+kk) - (a2/48) * Cb(9+kk) ) ;
    Cs = Sb(1+kk) - (a2/8)*( Sb(5+kk) - (a2/48) * Sb(9+kk) ) ;
    Sc = (a/2)*( Cb(3+kk) - (a2/24)*(Cb(7+kk)-(a2/80)*Cb(11+kk) ) ) ;
    Ss = (a/2)*( Sb(3+kk) - (a2/24)*(Sb(7+kk)-(a2/80)*Sb(11+kk) ) ) ;
    X(1+kk) = Cc - Ss ;
    Y(1+kk) = Sc + Cs ;
  end
end
%
%
%
function [X,Y] = GeneralizedFresnel_b_small( nk, a, b )

  X = zeros(nk,1) ;
  Y = zeros(nk,1) ;

  [Ca,Sa] = int_sincos2_tau_k( 5+nk, a, b ) ;
  
  % cos(b*tau) ~ 1-(1/2)*b^2*tau^2+(1/24)*b^4*tau^4 ;
  % sin(b*tau) ~ b*tau*(1-(1/6)*b^2*tau^2+(1/120)*b^4*tau^4) ;
  
  b2 = b*b ;
  for kk=0:nk-1
    Cc = Ca(1+kk) - (b2/2)*( Ca(3+kk) - (b2/12) * Ca(5+kk) ) ;
    Cs = b * ( Ca(2+kk) - (b2/6)*( Ca(4+kk) - (b2/20) * Ca(6+kk) ) ) ;
    Sc = Sa(1+kk) - (b2/2)*( Sa(3+kk) - (b2/12) * Sa(5+kk) ) ;
    Ss = b * ( Sa(2+kk) - (b2/6)*( Sa(4+kk) - (b2/20) * Sa(6+kk) ) ) ;
    X(1+kk) = Cc - Ss ;
    Y(1+kk) = Sc + Cs ;
  end
end
%
%
%
function [X,Y] = GeneralizedFresnel_ab_small( nk, a, b )

  X = zeros(nk,1) ;
  Y = zeros(nk,1) ;

  a2 = a*a ;
  a4 = a2*a2 ;
  b2 = b*b ;
  b4 = b2*b2 ;
  
  X(1) = a * b * ((-b2 / 0.2880e4 + 0.1e1 / 0.384e3) * a2 - b4 / 0.1920e4 + b2 / 0.72e2 - 0.1e1 / 0.8e1 - a4 / 0.46080e5) + a2 * (-b4 / 0.1728e4 + b2 / 0.112e3 - 0.1e1 / 0.40e2 - a4 / 0.599040e6) + (-a4 / 0.8448e4 - b4 / 0.5040e4 - 0.1e1 / 0.6e1) * b2 + b4 / 0.120e3 + 0.1e1 + a4 / 0.3456e4 ; 
  Y(1) = a * (a2 * (-b4 / 0.12672e5 - a4 / 0.9676800e7 + b2 / 0.864e3 - 0.1e1 / 0.336e3) + (-b4 / 0.12960e5 - a4 / 0.99840e5 - 0.1e1 / 0.20e2) * b2 + b4 / 0.336e3 + 0.1e1 / 0.6e1 + a4 / 0.42240e5) + b * (a2 * (-b4 / 0.9600e4 - a4 / 0.645120e6 + b2 / 0.384e3 - 0.1e1 / 0.48e2) + (-b4 / 0.40320e5 - a4 / 0.27648e5 - 0.1e1 / 0.24e2) * b2 + b4 / 0.720e3 + 0.1e1 / 0.2e1 + a4 / 0.3840e4) ; 

  if nk < 2 ; return ; end ;

  X(2) = a * b * ((-b2 / 0.3168e4 + 0.1e1 / 0.432e3) * a2 - b4 / 0.2160e4 + b2 / 0.84e2 - 0.1e1 / 0.10e2 - a4 / 0.49920e5) + a2 * (-b4 / 0.1920e4 + b2 / 0.128e3 - 0.1e1 / 0.48e2 - a4 / 0.645120e6) + (-a4 / 0.9216e4 - b4 / 0.5760e4 - 0.1e1 / 0.8e1) * b2 + b4 / 0.144e3 + 0.1e1 / 0.2e1 + a4 / 0.3840e4 ; 
  Y(2) = a * (a2 * (-b4 / 0.13824e5 - a4 / 0.10321920e8 + b2 / 0.960e3 - 0.1e1 / 0.384e3) + (-b4 / 0.14400e5 - a4 / 0.107520e6 - 0.1e1 / 0.24e2) * b2 + b4 / 0.384e3 + 0.1e1 / 0.8e1 + a4 / 0.46080e5) + b * (a2 * (-b4 / 0.10560e5 - a4 / 0.691200e6 + b2 / 0.432e3 - 0.1e1 / 0.56e2) + (-b4 / 0.45360e5 - a4 / 0.29952e5 - 0.1e1 / 0.30e2) * b2 + b4 / 0.840e3 + 0.1e1 / 0.3e1 + a4 / 0.4224e4);

  if nk < 3 ; return ; end ;

  X(3) = a * b * ((-b2 / 0.3456e4 + 0.1e1 / 0.480e3) * a2 - b4 / 0.2400e4 + b2 / 0.96e2 - 0.1e1 / 0.12e2 - a4 / 0.53760e5) + a2 * (-b4 / 0.2112e4 + b2 / 0.144e3 - 0.1e1 / 0.56e2 - a4 / 0.691200e6) + (-a4 / 0.9984e4 - b4 / 0.6480e4 - 0.1e1 / 0.10e2) * b2 + b4 / 0.168e3 + 0.1e1 / 0.3e1 + a4 / 0.4224e4;
  Y(3) = a * (a2 * (-b4 / 0.14976e5 - a4 / 0.10967040e8 + b2 / 0.1056e4 - 0.1e1 / 0.432e3) + (-b4 / 0.15840e5 - a4 / 0.115200e6 - 0.1e1 / 0.28e2) * b2 + b4 / 0.432e3 + 0.1e1 / 0.10e2 + a4 / 0.49920e5) + b * (a2 * (-b4 / 0.11520e5 - a4 / 0.737280e6 + b2 / 0.480e3 - 0.1e1 / 0.64e2) + (-b4 / 0.50400e5 - a4 / 0.32256e5 - 0.1e1 / 0.36e2) * b2 + b4 / 0.960e3 + 0.1e1 / 0.4e1 + a4 / 0.4608e4);

  if nk < 4 ; return ; end ;

  X(4) = a * b * ((-b2 / 0.3744e4 + 0.1e1 / 0.528e3) * a2 - b4 / 0.2640e4 + b2 / 0.108e3 - 0.1e1 / 0.14e2 - a4 / 0.57600e5) + a2 * (-b4 / 0.2304e4 + b2 / 0.160e3 - 0.1e1 / 0.64e2 - a4 / 0.737280e6) + (-a4 / 0.10752e5 - b4 / 0.7200e4 - 0.1e1 / 0.12e2) * b2 + b4 / 0.192e3 + 0.1e1 / 0.4e1 + a4 / 0.4608e4;
  Y(4) = a * (a2 * (-b4 / 0.16128e5 - a4 / 0.11612160e8 + b2 / 0.1152e4 - 0.1e1 / 0.480e3) + (-b4 / 0.17280e5 - a4 / 0.122880e6 - 0.1e1 / 0.32e2) * b2 + b4 / 0.480e3 + 0.1e1 / 0.12e2 + a4 / 0.53760e5) + b * (a2 * (-b4 / 0.12480e5 - a4 / 0.783360e6 + b2 / 0.528e3 - 0.1e1 / 0.72e2) + (-b4 / 0.55440e5 - a4 / 0.34560e5 - 0.1e1 / 0.42e2) * b2 + b4 / 0.1080e4 + 0.1e1 / 0.5e1 + a4 / 0.4992e4);

  if nk < 5 ; return ; end ;

  X(5) = a * b * ((-b2 / 0.4032e4 + 0.1e1 / 0.576e3) * a2 - b4 / 0.2880e4 + b2 / 0.120e3 - 0.1e1 / 0.16e2 - a4 / 0.61440e5) + a2 * (-b4 / 0.2496e4 + b2 / 0.176e3 - 0.1e1 / 0.72e2 - a4 / 0.783360e6) + (-a4 / 0.11520e5 - b4 / 0.7920e4 - 0.1e1 / 0.14e2) * b2 + b4 / 0.216e3 + 0.1e1 / 0.5e1 + a4 / 0.4992e4;
  Y(5) = a * (a2 * (-b4 / 0.17280e5 - a4 / 0.12257280e8 + b2 / 0.1248e4 - 0.1e1 / 0.528e3) + (-b4 / 0.18720e5 - a4 / 0.130560e6 - 0.1e1 / 0.36e2) * b2 + b4 / 0.528e3 + 0.1e1 / 0.14e2 + a4 / 0.57600e5) + b * (a2 * (-b4 / 0.13440e5 - a4 / 0.829440e6 + b2 / 0.576e3 - 0.1e1 / 0.80e2) + (-b4 / 0.60480e5 - a4 / 0.36864e5 - 0.1e1 / 0.48e2) * b2 + b4 / 0.1200e4 + 0.1e1 / 0.6e1 + a4 / 0.5376e4);

  if nk > 5
    error( 'nk > 5 ') ;
  end ;

end
%
%
%
function [X,Y] = GeneralizedFresnel_ab_large( nk, aa, b )
  a    = abs(aa) ;
  ras  = sqrt(pi / a);
  tmp  = sqrt(a*pi) ;
  [FC0,FS0] = FresnelCS(b/tmp);
  [FC1,FS1] = FresnelCS((a - b)/tmp);
  [FC2,FS2] = FresnelCS((a + b)/tmp);
  tmp  = b^2/(2*a) ;
  S    = sin(tmp);
  C    = cos(tmp);
  tmp  = a /2 + b ;
  Sapb = sin(tmp);
  Capb = cos(tmp);
  tmp  = a /2 - b ;
  Samb = sin(tmp);
  Camb = cos(tmp);
  X    = zeros(nk,1) ;
  Y    = zeros(nk,1) ;
  if aa > 0

    tmpC2 = FC2 - FC0 ;
    tmpS2 = FS2 - FS0 ;

    X(1) = ((tmpC2 * C + S * tmpS2) * ras);
    Y(1) = (tmpS2 * C - tmpC2 * S) * ras;

    if nk < 2 ; return ; end ;

    X(2) = ((-tmpC2 * C - S * tmpS2) * b * ras + Sapb) / a;
    Y(2) = ((-tmpS2 * C + tmpC2 * S) * b * ras - Capb + 1) / a;

    if nk < 3 ; return ; end ;

    X(3) = (ras * (-tmpS2 * C + tmpC2 * S) + Sapb) / a + ((tmpC2 * C + S * tmpS2) * ras * b ^ 2 - Sapb * b) / a ^ 2;
    Y(3) = ((tmpC2 * C + S * tmpS2) * ras - Capb) / a + ((tmpS2 * C - tmpC2 * S) * ras * b ^ 2 + b * (Capb - 1)) / a ^ 2;

    if nk < 4 ; return ; end ;

    X(4) = Sapb / a + (((3 * tmpS2 * C - 3 * tmpC2 * S) * ras - Sapb) * b + 2 * Capb - 2) / a ^ 2 + ((-tmpC2 * C - S * tmpS2) * ras * b ^ 3 + Sapb * b ^ 2) / a ^ 3;
    Y(4) = -Capb / a + (((-3 * tmpC2 * C - 3 * S * tmpS2) * ras + Capb) * b + 2 * Sapb) / a ^ 2 + ((-tmpS2 * C + tmpC2 * S) * ras * b ^ 3 + (-Capb + 1) * b ^ 2) / a ^ 3;

    if nk < 5 ; return ; end ;

    X(5) = Sapb / a + (-Sapb * b + (-3 * tmpC2 * C - 3 * S * tmpS2) * ras + 3 * Capb) / a ^ 2 + (((-6 * tmpS2 * C + 6 * tmpC2 * S) * ras + Sapb) * b ^ 2 + (-5 * Capb + 5) * b) / a ^ 3 + ((tmpC2 * C + S * tmpS2) * ras * b ^ 4 - Sapb * b ^ 3) / a ^ 4;
    Y(5) = -Capb / a + (Capb * b + (-3 * tmpS2 * C + 3 * tmpC2 * S) * ras + 3 * Sapb) / a ^ 2 + (((6 * tmpC2 * C + 6 * S * tmpS2) * ras - Capb) * b ^ 2 - 5 * Sapb * b) / a ^ 3 + ((tmpS2 * C - tmpC2 * S) * ras * b ^ 4 + (Capb - 1) * b ^ 3) / a ^ 4;

  else

    tmpC1 = FC1 + FC0 ;
    tmpS1 = FS1 + FS0 ;

    X(1) = ras * (tmpC1 * C + tmpS1 * S);
    Y(1) = ((-tmpS1 * C + tmpC1 * S) * ras);

    if nk < 2 ; return ; end ;

    X(2) = (b * (tmpC1 * C + tmpS1 * S) * ras + Samb) / a;
    Y(2) = ((-tmpS1 * C + tmpC1 * S) * ras * b + Camb - 1) / a;

    if nk < 3 ; return ; end ;

    X(3) = ((-tmpS1 * C + tmpC1 * S) * ras + Samb) / a + (b ^ 2 * (tmpC1 * C + tmpS1 * S) * ras + Samb * b) / a ^ 2;
    Y(3) = ((-tmpC1 * C - tmpS1 * S) * ras + Camb) / a + ((-tmpS1 * C + tmpC1 * S) * ras * b ^ 2 + b * (Camb - 1)) / a ^ 2;

    if nk < 4 ; return ; end ;

    X(4) = Samb / a + (((-3 * tmpS1 * C + 3 * tmpC1 * S) * ras + Samb) * b + 2 * Camb - 2) / a ^ 2 + (ras * (tmpC1 * C + tmpS1 * S) * b ^ 3 + Samb * b ^ 2) / a ^ 3;
    Y(4) = Camb / a + (((-3 * tmpC1 * C - 3 * tmpS1 * S) * ras + Camb) * b - 2 * Samb) / a ^ 2 + ((-tmpS1 * C + tmpC1 * S) * ras * b ^ 3 + (Camb - 1) * b ^ 2) / a ^ 3;

    if nk < 5 ; return ; end ;

    X(5) = Samb / a + (Samb * b + (-3 * tmpC1 * C - 3 * tmpS1 * S) * ras + 3 * Camb) / a ^ 2 + (((-6 * tmpS1 * C + 6 * tmpC1 * S) * ras + Samb) * b ^ 2 + (5 * Camb - 5) * b) / a ^ 3 + (ras * (tmpC1 * C + tmpS1 * S) * b ^ 4 + Samb * b ^ 3) / a ^ 4;
    Y(5) = Camb / a + (Camb * b + (3 * tmpS1 * C - 3 * tmpC1 * S) * ras - 3 * Samb) / a ^ 2 + (((-6 * tmpC1 * C - 6 * tmpS1 * S) * ras + Camb) * b ^ 2 - 5 * Samb * b) / a ^ 3 + ((-tmpS1 * C + tmpC1 * S) * ras * b ^ 4 + (Camb - 1) * b ^ 3) / a ^ 4;

  end
end
