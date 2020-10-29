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
y0     = 200;
theta0 = pi/12;
kappa0 = 1/6179.95194205085;
L0     = 201.704876743562; % 220.840282897566;

C0 = CircleArc( x0, y0, theta0, kappa0, L0 );

x1     = 400;
y1     = 200;
theta1 = -pi/4;
kappa1 = -1/199.418998889075;
L1     = -220.840282897566; %201.704876743562;
C1     = CircleArc( x1, y1, theta1, kappa1, L1 );


hold off
C0.plot(100,{'Color','red','Linewidth',3});
hold on
C1.plot(100,{'Color','blue','Linewidth',3});
%C2.plot(100,{'Color','black','Linewidth',2});
axis equal

L1 = Biarc(x0,y0,theta0,x1,y1,theta1);
L1.info();
L1.plot();
