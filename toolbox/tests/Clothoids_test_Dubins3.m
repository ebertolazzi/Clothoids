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

DB1 = Dubins();
DB2 = Dubins();

k_max = 1;
  
x0     = -0.011;
y0     = 0.2;
x3     = x0+1;
y3     = y0;
theta0 = -pi/2;
theta3 = pi/2;
DB1.build( x0, y0, theta0, x3, y3, theta3, k_max );
DB1.plot();

hold on;

x0     = 1;
y0     = 0;
x3     = x0+1;
y3     = y0;
theta0 = pi/2;
theta3 = pi/2;
DB2.build( x0, y0, theta0, x3, y3, theta3, k_max );
DB2.plot();

axis equal;

[s1,s2] = DB1.intersect( DB2 );

XY1 = DB1.eval( s1 );
XY2 = DB2.eval( s2 );

plot( XY1(1,:), XY1(2,:), 'ro', 'MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor',[0.5,0.5,0.5]);
plot( XY2(1,:), XY2(2,:), 'b*', 'MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor',[0.9,0.9,0.5]);


