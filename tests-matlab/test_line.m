addpath('../matlab');

% check constructors
x0     = 0 ;
y0     = 2 ;
x1     = 4 ;
y1     = 3 ;
theta0 = 0 ;
L      = 10 ;
L1 = LineSegment( x0, y0, theta0, L );
L2 = LineSegment( x0, y0, theta0+pi/4, 0, L );
L3 = LineSegment( [x0, y0+1], [x1, y1] );

L1.plot() ;
L2.plot('red',3) ;
L3.plot('blue',3) ;

axis equal;
