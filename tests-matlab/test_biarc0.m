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

addpath('../matlab');
close all;

% check constructors
x0     = 0;
y0     = 2;
x1     = 3;
y1     = 3;
x2     = 5;
y2     = 2;
theta0 = 0;
theta1 = pi;
L      = 10;
k0     = 1/3;
dk     = 0;

L1 = Biarc( x0, y0, theta0, x1, y1, theta1 );
L2 = Biarc();
L2.build( x0, y0, theta0, x1, y1, theta1 );
L3 = Biarc();
L3.build( x0, y0, theta0, x2, y2, theta1 );

%
npts = 1000;
fmt1 = {'Color','red','LineWidth',3};
fmt2 = {'Color','black','LineWidth',3};
fmt3 = {'Color','green','LineWidth',3};
L1.plot(npts,fmt1);
hold on;
L2.plot(npts,fmt1);
L3.plot(npts,fmt2);

L1.translate(1,1);
L1.plot();

L3.rotate(pi/3,0,0);
L3.plot(npts,fmt1,fmt3);

L3.rotate(pi/3,0,0);
L3.plot(npts,fmt1,fmt3);

L3.rotate(pi/3,0,0);
L3.plot(npts,fmt1,fmt3);

L3.translate(3,3);
L3.plot(npts,fmt1,fmt3);

L3.eval(1)
L3.eval_D(1)
L3.eval_DD(1)
L3.eval_DDD(1)
L3.xBegin()
L3.yBegin()
L3.thetaBegin()
L3.length()

%L3.distance(1,4)

axis equal;
%