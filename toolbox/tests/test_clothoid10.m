%=========================================================================%
%                                                                         %
%  Autors: Enrico Bertolazzi                                              %
%          Department of Industrial Engineering                           %
%          University of Trento                                           %
%          enrico.bertolazzi@unitn.it                                     %
%          m.fregox@gmail.com                                             %
%                                                                         %
%=========================================================================%
% Driver test program to check Clothoids lib                              %
%=========================================================================%

close all;

% check constructors
x0     = -5;
y0     = 10;
theta0 = 0;
kappa0 = -0.6;
dk     = 0.1;
L      = 15;

L1 = ClothoidCurve( x0, y0, theta0, kappa0, dk, L );
L2 = ClothoidCurve();
L2.copy(L1);
L2.trim(1.23,10);

[x1,y1,theta1,kappa1] = L1.evaluate( 1.23+0.1 );
fprintf('x = %g y = %g theta = %g kappa = %g\n', x1, y1, theta1, kappa1 );
[x2,y2,theta2,kappa2] = L2.evaluate( 0+0.1 );
fprintf('x = %g y = %g theta = %g kappa = %g\n', x2, y2, theta2, kappa2 );

fprintf('dx = %g dy = %g dtheta = %g dkappa = %g\n', x1-x2, y1-y2, theta1-theta2, kappa1-kappa2 );
