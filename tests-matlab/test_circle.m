addpath('../matlab');

C1 = CircleArc( 0,0,0,0.1,12);
C2 = CircleArc( 0,0,0,-0.2,12);

C1.plot() ;
C2.plot('red',3) ;
axis equal;
