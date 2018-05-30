addpath('../matlab');
close all;

% check constructors
x0     = 0 ;
y0     = 2 ;
x1     = 3 ;
y1     = 3 ;
x2     = 5 ;
y2     = 2 ;
theta0 = 0 ;
theta1 = pi;
L      = 10 ;
k0     = 1/3 ;
dk     = 0;

L1 = Biarc( x0, y0, theta0, x1, y1, theta1 );
L2 = Biarc() ;
L2.build( x0, y0, theta0, x1, y1, theta1 );
L3 = Biarc() ;
L3.build( x0, y0, theta0, x2, y2, theta1 );

%
step = 0.01;
L1.plot(step,'Color','red','LineWidth',3) ;
hold on ;
L2.plot(step,'Color','red','LineWidth',3) ;
L3.plot(step,'Color','black','LineWidth',3) ;

L1.translate(1,1) ;
L1.plot() ;

L3.rotate(pi/3,0,0) ;
L3.plot(step,'Color','green','LineWidth',3) ;

L3.rotate(pi/3,0,0) ;
L3.plot(step,'Color','green','LineWidth',3) ;

L3.rotate(pi/3,0,0) ;
L3.plot(step,'Color','green','LineWidth',3) ;

L3.translate(3,3) ;
L3.plot(step,'Color','blue','LineWidth',3) ;

L3.eval(1)
L3.eval_D(1)
L3.eval_DD(1)
L3.eval_DDD(1)
L3.xBegin0()
L3.yBegin0()
L3.thetaBegin0()
L3.length()

%L3.distance(1,4)

axis equal;
%