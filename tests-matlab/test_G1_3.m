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

arange = pi*0.99;
brange = pi*0.5;
rangle = pi/2; % rotate curves by rangle

thetaV = [-brange:brange/5:brange];
phiV   = [-arange:arange/6:arange];

x0     = 0; % initial point with initial angle
y0     = 0;
theta0 = rangle;

fprintf('Testing G1 Clothoid interpolation on multiple points\n');
fprintf('fixed initial point (%g,%g) initial angle = %g\n', x0, y0, theta0);

hFig = figure(1);
set(hFig, 'Position', [1 1 1000 800]);

color = { [0.5 0.4 1], [0.9 0 0], [1 0.5 0], [0.6 0.8 0.2] , [0.9 0.9 0.2] };
kk = 1;
hold on;
C = ClothoidCurve();
for theta=thetaV
  for phi=phiV
    x1     = cos(phi+rangle);
    y1     = sin(phi+rangle);
    theta1 = phi+rangle+theta;
    fprintf('final point (%g,%g) initial angle = %g\n', x1, y1, theta1);
    C.build_G1( x0, y0, theta0, x1, y1, theta1 );
    C.plot( 1000, '-','Color', color{kk} ); % plot computed curve
  end
  kk = kk+1;
  if kk > 5
    kk = 1;
  end
end
hold off;
axis equal;
