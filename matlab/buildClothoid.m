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
function [S,varargout] = buildClothoid( x0, y0, theta0, x1, y1, theta1, varargin )

  dx  = x1 - x0 ;
  dy  = y1 - y0 ;
  r   = sqrt( dx^2 + dy^2 ) ;
  phi = atan2( dy, dx ) ;

  phi0  = normalizeAngle(theta0 - phi) ;
  phi1  = normalizeAngle(theta1 - phi) ;
  delta = phi1 - phi0 ;

  % initial point
  Aguess = guessA( phi0, phi1 ) ;

  % Newton iteration
  [A,iter] = findA( Aguess, delta, phi0, 1e-12 ) ;
  if nargout == 2
    varargout{2} = iter ;
  end

  % final operation
  [h,g] = GeneralizedFresnelCS( 1, 2*A, delta-A, phi0 ) ;
  L = r/h ;

  if L > 0
    k  = (delta - A)/L ;
    dk = 2*A/L^2 ;
  else
    error('negative length') ;
  end

  S.x     = x0 ;
  S.y     = y0 ;
  S.theta = theta0 ;
  S.k     = k ;
  S.dk    = dk ;
  S.L     = L ;
  
  if nargin == 7 && varargin{1}

    [X,Y] = GeneralizedFresnelCS( 3, 2*A, delta-A, theta0 ) ;
    
    if true
      alpha = X(1)*X(2) + Y(1)*Y(2) ;
      beta  = X(1)*X(3) + Y(1)*Y(3) ;
      gamma = X(1)^2+Y(1)^2 ;
      tx    = X(2)-X(3) ;
      ty    = Y(2)-Y(3) ;
      txy   = L*(X(2)*Y(3)-X(3)*Y(2)) ;
      omega = L*(Y(1)*tx-X(1)*ty) - txy ;
      delta = X(1)*tx + Y(1)*ty ;

      L_1  = omega/delta ; % L_0
      L_2  = txy/delta ; % L_1

      delta = delta * L ;
      k_1  = (beta-gamma-k*omega)/delta ; % k_0
      k_2  = -(beta+k*txy)/delta ; % k_1

      delta = delta * L/2 ;
      dk_1 = (gamma-alpha-dk*omega*L)/delta ; % dk_0    
      dk_2 = (alpha-dk*txy*L)/delta ; % dk_1
    else
      dkL = dk*L ;
      M = [ X(1) - L*(dkL*Y(3)+k*Y(2)), -L*Y(2), -L*Y(3) ; ...
            Y(1) + L*(dkL*X(3)+k*X(2)),  L*X(2),  L*X(3) ; ...
            dkL+k,                            1,       1 ] ;
      tmp = M\[L*Y(1);-L*X(1);-1] ;
      L_1  = tmp(1) ;
      k_1  = tmp(2)/L ;
      dk_1 = 2*tmp(3)/L^2 ;
      tmp = M\[0;0;1] ;
      L_2  = tmp(1) ;
      k_2  = tmp(2)/L ;
      dk_2 = 2*tmp(3)/L^2 ;
    end
    
    S.k_1  = k_1 ; % k_0
    S.dk_1 = dk_1 ; % dk_0
    S.L_1  = L_1 ; % L_0
    
    S.k_2  = k_2 ; % k_1
    S.dk_2 = dk_2 ; % dk_1
    S.L_2  = L_2 ; % L_1
  end
  
end

%=========================================================================%
%  normalizeAngle:  normalize angle in the range [-pi,pi]                 %
%=========================================================================%
function phi = normalizeAngle( phi_in )
  phi = phi_in ;
  while ( phi > pi )
    phi = phi - 2*pi ;
  end
  while ( phi < -pi )
    phi = phi + 2*pi ;
  end
end

%=========================================================================%
%  findA:  Find a zero of function g(A) defined as                        %
%  g(A) = \int_0^1 \sin( A*t^2+(delta-A)*t+phi0 ) dt                      %
%                                                                         %
%  USAGE:  A = findA( Aguess, delta, phi0, tol );                         %
%                                                                         %
%  Given an initial guess Aguess find the closest zero of equation g(A)   %
%                                                                         %
%  On input:                                                              %
%    Aguess      = initial guess.                                         %
%    delta, phi0 = Angles used in the clothoid fitting problem.           %
%    tol         = Tolerance for stopping criterium of Newton iteration.  %
%                                                                         %
%  On output:                                                             %
%    A           = the zero of function g(A) closest to Aguess.           %
%    iter        = iteration performed                                    %
%                                                                         %
%=========================================================================%
function [A,iter] = findA( Aguess, delta, phi0, tol )
  A = Aguess ;
  for iter=1:100
    [intC,intS] = GeneralizedFresnelCS( 3, 2*A, delta-A, phi0 ) ;    
    f  = intS(1) ;
    df = intC(3)-intC(2) ;
    A  = A - f/df ;
    if abs(f) < tol
      break ;
    end
  end
  if abs(f) > tol*10
    fprintf( 1, 'Newton iteration fails, f = %g\n', f ) ;
    fprintf( 1, 'Aguess = %g, A = %g, delta = %g , phi0 = %g\n', Aguess, A, delta, phi0 ) ;
  end
end

%=========================================================================%
%  guessA:  Find guess for zeros of function g(A)                         %
%                                                                         %
%  USAGE:  A = guessA( phi0, phi1 );                                      %
%                                                                         %
%  On input:                                                              %
%       phi0, phi1 = Angles used in the clothoid fitting problem.         %
%                                                                         %
%  On output:                                                             %
%       A = an approximate zero of function g(A).                         %
%                                                                         %
%=========================================================================%
function A = guessA( phi0, phi1 )
  CF = [ 2.989696028701907, ...
         0.716228953608281, ...
        -0.458969738821509, ...
        -0.502821153340377, ...
         0.261062141752652, ...
        -0.045854475238709 ] ;
  X  = phi0/pi ;
  Y  = phi1/pi ;
  xy = X*Y ;
  A  = (phi0+phi1) * ( CF(1) + xy * ( CF(2) + xy * CF(3) ) + ...
                       (CF(4) + xy * CF(5)) * (X^2+Y^2) + ...
                        CF(6) * (X^4+Y^4) ) ;
end
