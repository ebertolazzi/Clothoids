addpath('../matlab');

%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi                                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
% Driver test program to check bounding box on clothoid                       %
%=============================================================================%

close all ;

x0     = 0 ;
y0     = 0 ;
theta0 = pi*0.7 ;
k0     = -0.1 ;
dk     = 0.05 ;
L      = 20 ;
C1 = ClothoidCurve( x0, y0, theta0, k0, dk, L );

x0     = -10 ;
y0     = 0 ;
theta0 = pi/4 ;
k0     = 0.1 ;
dk     = -0.05 ;
L      = 20 ;
C2 = ClothoidCurve( x0, y0, theta0, k0, dk, L );

max_angle = pi/2 ;
max_size  = 0.5 ;

step = 0.01 ;

C1.plot(step,'Color','red','LineWidth',3) ;
hold on ;
C2.plot(step,'Color','black','LineWidth',3) ;

[s1,s2] = intersectClothoid( C1, C2 ) ;

s1
s2

XY1 = C1.eval( s1 )
XY2 = C2.eval( s2 )

plot( XY1(1,:), XY1(2,:), 'ob', 'LineWidth', 2 ) ;
plot( XY2(1,:), XY2(2,:), 'oc', 'LineWidth', 2 ) ;

axis equal ;
