%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi and Marco Frego                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
% Driver test program to check clothoid computation                           %
%=============================================================================%

close all ;

npts = 400 ;

% initial point with angle direction
x0     = 0 ;
y0     = 0.1 ;
theta0 = 1.571 ;

% final point with angle direction
x1     = 0.085 ;
y1     = 0.035 ;
theta1 = 2.356 ;

fprintf('Testing G1 Clothoid interpolation\n') ;
fprintf('initial point (%g,%g) initial angle = %g\n', x0, y0, theta0) ;
fprintf('final point (%g,%g) final angle = %g\n', x1, y1, theta1) ;

% compute clothoid parameters
S = buildClothoid( x0, y0, theta0, x1, y1, theta1 ) ;

fprintf('Computed parameters: k = %g, k'' = %g, L = %g\n', S.k, S.dk, S.L ) ;

% compute points on clothoid
XY = pointsOnClothoid( S, 0:S.L/npts:S.L ) ;

plot( XY(1,:), XY(2,:), '-r', 'LineWidth', 3 ) ;
axis equal
