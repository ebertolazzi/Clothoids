addpath('../matlab');

close all ;
clear all ;

% check constructors
x0     = [-5,0,0,2.5] ;
y0     = [10,0,1,2] ;
theta0 = 0 ;
kappa0 = [-0.6, 0.2, 0.2, 0.5] ;
L      = [5,30,100,30] ;


for kk=1:4
  subplot(2,2,kk) ;
  L1 = CircleArc( x0(kk), y0(kk), theta0, kappa0(kk), L(kk) );
  %
  L1.plot() ;
  hold on;

  x     = -10:0.05:10 ;
  y     = -5:0.05:15 ;
  [X,Y] = meshgrid(x,y);

  tic
  Z = L1.distance(X,Y);
  toc

  contour(X,Y,Z,100)
  %surf(X,Y,Z)
  axis equal;
  %
  L1.delete() ;
end

