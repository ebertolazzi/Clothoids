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

addpath('../G1fitting');

close all ;

c1.x0        = 0 ;
c1.y0        = 0 ;
c1.theta0    = pi*0.7 ;
c1.kappa     = -0.1 ;
c1.dkappa    = 0.15 ;
c1.L         = 20 ;

c2.x0        = -10 ;
c2.y0        = 0 ;
c2.theta0    = pi/4 ;
c2.kappa     = 0.1 ;
c2.dkappa    = -0.05 ;
c2.L         = 20 ;

max_angle = pi/2 ;
max_size  = 0.5 ;

npts = 10000 ;

XY1 = pointsOnClothoid( c1.x0, c1.y0, c1.theta0, c1.kappa, c1.dkappa, c1.L, npts ) ;
XY2 = pointsOnClothoid( c2.x0, c2.y0, c2.theta0, c2.kappa, c2.dkappa, c2.L, npts ) ;

offs = 0 ;

subplot(3,1,1) ;
hold on ;
TT1 = bbClothoid( c1.x0, c1.y0, c1.theta0, c1.kappa, c1.dkappa, c1.L, max_angle, max_size, offs ) ;
TT2 = bbClothoid( c2.x0, c2.y0, c2.theta0, c2.kappa, c2.dkappa, c2.L, max_angle, max_size, offs ) ;
size(TT1)
size(TT2)

for i=1:size(TT1,2)
  fill( TT1(1:2:end,i), TT1(2:2:end,i), 'red') ;
end

for i=1:size(TT2,2)
  fill( TT2(1:2:end,i), TT2(2:2:end,i), 'yellow') ;
end

plot( XY1(1,:), XY1(2,:), '-b', 'LineWidth', 0.1 ) ;
plot( XY2(1,:), XY2(2,:), '-k', 'LineWidth', 0.1 ) ;

%AXY1 = pointsOnClothoid( c1.x0, c1.y0, c1.theta0, c1.kappa, c1.dkappa, c1.L/2 ) ;
%AXY2 = pointsOnClothoid( c2.x0, c2.y0, c2.theta0, c2.kappa, c2.dkappa, c2.L/2 ) ;
%plot( AXY1(1,:), AXY1(2,:), 'ob', 'LineWidth', 0.1 ) ;
%plot( AXY2(1,:), AXY2(2,:), 'ok', 'LineWidth', 0.1 ) ;


axis equal ;

subplot(3,1,2) ;
hold on ;

plot( XY1(1,:), XY1(2,:), '-b', 'LineWidth', 1 ) ;
plot( XY2(1,:), XY2(2,:), '-k', 'LineWidth', 1 ) ;

for i=1:size(TT1,2)
  p0 = [ TT1(1,i), TT1(2,i) ] ;
  p1 = [ TT1(3,i), TT1(4,i) ] ;
  p2 = [ TT1(5,i), TT1(6,i) ] ;
  for j=1:size(TT2,2)
    q0 = [ TT2(1,j), TT2(2,j) ] ;
    q1 = [ TT2(3,j), TT2(4,j) ] ;
    q2 = [ TT2(5,j), TT2(6,j) ] ;
    ii = TriTriOverlap( p0, p1, p2, q0, q1, q2 ) ; 
    if ii ~= 0
      fill( TT1(1:2:end,i), TT1(2:2:end,i), 'red') ;
      fill( TT2(1:2:end,j), TT2(2:2:end,j), 'yellow') ;
    end
  end
end

axis equal ;


subplot(3,1,3) ;
hold on ;

plot( XY1(1,:), XY1(2,:), '-b', 'LineWidth', 1 ) ;
plot( XY2(1,:), XY2(2,:), '-k', 'LineWidth', 1 ) ;


[s1,s2] = intersectClothoid( c1, c2 ) ;

XY1 = pointsOnClothoid( c1.x0, c1.y0, c1.theta0, c1.kappa, c1.dkappa, s1 ) ;
XY2 = pointsOnClothoid( c2.x0, c2.y0, c2.theta0, c2.kappa, c2.dkappa, s2 ) ;

plot( XY1(1,:), XY1(2,:), 'ob', 'LineWidth', 2 ) ;
plot( XY2(1,:), XY2(2,:), 'ok', 'LineWidth', 2 ) ;


%for offs=[-0.5,0,0.5]
%  TT = bbClothoid( x0, y0, theta0, kappa, dkappa, L, max_angle, max_size, offs ) ;
  %for i=1:size(TT,2)
    %fill( TT(1:2:end,i), TT(2:2:end,i), 'red') ;
    %end
  %end
%plot( XY(1,:), XY(2,:), '-b', 'LineWidth', 2 ) ;

axis equal ;
