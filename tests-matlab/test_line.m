addpath('../matlab');

% check constructors
x0     = 0 ;
y0     = 2 ;
x1     = 4 ;
y1     = 3.5 ;
theta0 = 0 ;
L      = 10 ;
L1 = LineSegment( x0, y0, theta0, L );
L2 = LineSegment( x0, y0, theta0+pi/4, L );
L3 = LineSegment( [x0, y0+1], [x1, y1] );

L1.plot() ;
L2.plot('red',3) ;
L3.plot('blue',3) ;

L1.translate(1,1) ;
L1.plot() ;

L2.changeOrigin(1,2) ;
L2.plot('red',3) ;


L3.trim(2,3) ;
L3.plot('blue',3) ;

L3.rotate(pi/4,3,3) ;
L3.plot('green',3) ;

L3.rotate(pi/4,3,3) ;
L3.plot('black',3) ;

L3.rotate(pi/4,3,3) ;
L3.plot('yellow',3) ;

L3.rotate(pi/4,3,3) ;
L3.plot('blue',3) ;


L3.translate(3,3) ;
L3.plot('blue',3) ;

L3.eval(1)
L3.eval_D(1)
L3.eval_DD(1)
L3.eval_DDD(1)
L3.xBegin()
L3.yBegin()
L3.theta()
L3.length()

L1.delete()
L2.delete()
L3.delete()

axis equal;
