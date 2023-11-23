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

DB = Dubins();
hold on;

%tiledlayout(2,2);%, 'Padding', 'none', 'TileSpacing', 'compact');
DY = 4;

offs = [0, 0, 0 ];

x0     = 0;
y0     = 0;
theta0 = pi/2;
x3     = 6;
y3     = 0;
theta3 = pi/2;
k_max  = 1;
r_min  = 1/k_max;

dubConnObj = dubinsConnection('MinTurningRadius',r_min);

DB.build( x0, y0, theta0, x3, y3, theta3, k_max );
DB.plot();

startPose = [x0 y0 theta0]+offs;
goalPose  = [x3 y3 theta3]+offs;
[pathSegObj, pathCosts] = connect(dubConnObj,startPose,goalPose);
show(pathSegObj{1})

axis equal;
