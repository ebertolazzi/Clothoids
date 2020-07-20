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

% check constructors
x0     = 0;
y0     = 2;
x1     = 3;
y1     = 3;
x2     = 5;
y2     = 2;
theta0 = 0;
L      = 10;
k0     = 1/3;
L1 = CircleArc( x0, y0, theta0, k0, L );
L2 = CircleArc( x0, y0, theta0+pi/4, 1.5*k0, L );
L3 = CircleArc();
L3.build_3P(x0, y0, x1, y1, x2, y2 );
%
L1.plot();
hold on;
fmt1 = {'Color','red','Linewidth',3};
fmt2 = {'Color','black','Linewidth',3};
fmt3 = {'Color','blue','Linewidth',3};
fmt4 = {'Color','green','Linewidth',3};
L2.plot(100,fmt1);
L3.plot(100,fmt2);

L1.translate(1,1);
L1.plot();

L2.changeCurvilinearOrigin(4,10);
L2.plot(100,fmt1);

L2.changeOrigin(1,2);
L2.plot(100,fmt1);

L3.trim(2,7);
L3.plot(100,fmt3);

L3.rotate(pi/3,0,0);
L3.plot(100,fmt3);

L3.rotate(pi/3,0,0);
L3.plot(100,fmt3);

L3.rotate(pi/3,0,0);
L3.plot(100,fmt3);

L3.translate(3,3);
L3.plot(100,fmt4);

L3.eval(1)
L3.eval_D(1)
L3.eval_DD(1)
L3.eval_DDD(1)
L3.xBegin()
L3.yBegin()
L3.thetaBegin()
L3.length()

L3.distance(1,4)

axis equal;
%