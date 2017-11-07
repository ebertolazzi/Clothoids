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

x0        = 0 ;
y0        = 0 ;
theta0    = 0 ;
kappa     = -1 ;
dkappa    = 0.1 ;
L         = 20 ;
max_angle = pi/2 ;
max_size  = 0.25 ;

XY = pointsOnClothoid( x0, y0, theta0, kappa, dkappa, 0:L/100:L ) ;
hold off ;
plot( XY(1,:), XY(2,:), '-b', 'LineWidth', 2 ) ;
hold on ;

for offs=[-0.5,0,0.5]
  TT = bbClothoid( x0, y0, theta0, kappa, dkappa, L, max_angle, max_size, offs ) ;
  for i=1:size(TT,2)
    fill( TT(1:2:end,i), TT(2:2:end,i), 'red') ;
  end
end
plot( XY(1,:), XY(2,:), '-b', 'LineWidth', 2 ) ;

axis equal ;
