addpath('../matlab');

% check constructors
x0     = 0 ;
y0     = 2 ;
x1     = 4 ;
y1     = 3 ;
x2     = 5 ;
y2     = 2 ;
theta0 = 0 ;
L      = 10 ;
k0     = 1/3 ;
L1 = CircleArc( x0, y0, theta0, k0, L );
%
L1.plot() ;

x = -10:0.05:10 ;
y = -5:0.05:15 ;
[X,Y] = meshgrid(x,y);

tic
for i=1:size(X,1)
  for j=1:size(X,2)
    Z(i,j) = L1.distance(X(i,j),Y(i,j));
  end
end

toc

contour(X,Y,Z,100)
%surf(X,Y,Z)
axis equal;
%
L1.delete() ;
