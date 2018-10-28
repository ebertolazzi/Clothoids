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

x0        = 0;
y0        = 0;
theta0    = 0;
kappa     = -1;
dkappa    = 0.1;
L         = 20;

CLOT = ClothoidCurve( x0, y0, theta0, kappa, dkappa, L );

npts = 1000;
CLOT.plot(npts,'Color','blue','LineWidth',2);
hold on;

max_angle = pi/2;
max_size  = 0.25;

for offs=[-0.5,0,0.5]
  TT = CLOT.bbox( max_angle, max_size, offs );
  for i=1:size(TT,2)
    fill( TT(1:2:end,i), TT(2:2:end,i), 'red');
  end
end
%plot( XY(1,:), XY(2,:), '-b', 'LineWidth', 2 );

axis equal;
