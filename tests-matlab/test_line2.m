addpath('../matlab');

% check constructors
x0     = 0 ;
y0     = 2 ;
x1     = 4 ;
y1     = 3 ;
x2     = 5 ;
y2     = 2 ;
theta0 = 0 ;
L      = 5 ;
L1 = LineSegment( x0, y0, theta0, L );
%
L1.plot() ;

x = -10:0.05:10 ;
y = -5:0.05:15 ;
[X,Y] = meshgrid(x,y);

tic
Z = L1.distance(X,Y);
toc

contour(X,Y,Z,100)
%surf(X,Y,Z)
axis equal;
%
L1.delete() ;
