addpath('../matlab');

% check constructors
x0     = 0 ;
y0     = 2 ;
x1     = 4 ;
y1     = 3 ;
x2     = 5 ;
y2     = 2 ;
theta0 = 0 ;
kappa0 = -0.6 ;
dk     = 0.2 ;
L      = 10 ;
L1 = ClothoidCurve( x0, y0, theta0, kappa0, dk, L );
%
L1.plot(L/1000,'Color','green','LineWidth',3) ;
hold on;

x     = -10:0.05:10 ;
y     = -5:0.05:15 ;
[X,Y] = meshgrid(x,y);

%W     = L1.distance(-10,1.55);

tic
Z = L1.distance(X,Y);
toc

contour(X,Y,Z,100)
%surf(X,Y,Z)
axis equal;
%
L1.delete() ;
