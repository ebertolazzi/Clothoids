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

clear all ;
close all ;

% G2 fitting data
x0     = -0.5 ;
y0     = -1  ;
theta0 = pi/2 ;
kappa0 = -0.1 ; %0.2 + 1e-10 ;

x1     = 0.5 ;
y1     = 1 ;
theta1 = pi/2 ;
kappa1 = 0; %0.2 + 0 ;

figure('Position',[1,1,800,500]);
subplot(2,2,[1,3]);

% compute G2 clothoid
[ S0, S1, SM, SG, iter ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1);
iter
if iter > 0 
  draw3curve( S0, S1, SM, false );
end
axis equal;

[X,Y]    = pointsOnClothoid( SG, 0:SG.L/100:SG.L );
plot(X,Y,'--k','Linewidth',2);
hold on

s0 = S0.L ;
s1 = S1.L ;
[X,Y] = pointsOnCircle( x0, y0, theta0, kappa0, (0:s0/100:s0) ) ; 
plot(X,Y,'-g','Linewidth',3);
hold on ;
xx0 = X(end) ;
yy0 = Y(end) ;
tth0 = theta0 + s0 * kappa0 ;

[X,Y] = pointsOnCircle( x1, y1, theta1, kappa1, -(0:s1/100:s1) ) ; 
plot(X,Y,'-g','Linewidth',3);
xx1  = X(end) ;
yy1  = Y(end) ;
tth1 = theta1 - s1 * kappa1 ;

subplot(2,2,2);
draw3curvature( S0, S1, SM, false );

subplot(2,2,4);
draw3angle( S0, S1, SM, false );
