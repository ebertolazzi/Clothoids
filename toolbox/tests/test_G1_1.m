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

% initial point with angle direction
x0     = 0.5;
y0     = 0;
theta0 = 1;

% final point with angle direction
x1     = 0.5;
y1     = 0.5;
theta1 = 1;

fprintf('Testing G1 Clothoid interpolation\n');
fprintf('initial point (%g,%g) initial angle = %g\n', x0, y0, theta0);
fprintf('final point (%g,%g) final angle = %g\n', x1, y1, theta1);

% compute clothoid parameters
C = ClothoidCurve();
C.build_G1( x0, y0, theta0, x1, y1, theta1 );

fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', C.kappaBegin(), C.dkappa(), C.length() );

C.plot();
axis equal
