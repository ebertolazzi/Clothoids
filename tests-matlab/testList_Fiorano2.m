addpath('../matlab');

%=========================================================================%
%                                                                         %
%  Autor: Enrico Bertolazzi                                               %
%         Department of Industrial Engineering                            %
%         University of Trento                                            %
%         enrico.bertolazzi@unitn.it                                      %
%                                                                         %
%=========================================================================%
% Driver test program to check Clothoids lib                              %
%=========================================================================%

close all;

SL = ClothoidList();

data   = importdata('data_logged_giro_Veloce.txt');
x      = data.data(:,3);
y      = data.data(:,4);

data   = importdata('fiorano-circuit-3D-kerbs-0.5m.txt');
s      = data.data(:,1);
kappa  = data.data(:,2);
x0     = 0;
y0     = 0;
theta0 = 0*pi;
ok     = SL.build( x0, y0, theta0, s, kappa);
SL.plot();
hold on;
plot( x, y, 'ob', 'LineWidth', 2 );

dst = SL.distance( x, y );
[ xx, yy, s, t, iflag, dst ] = SL.closestPoint( x, y );

%[ ss, tt, ipos ] = SL.find_coord1( x, y );
[ ss, tt ] = SL.find_coord( x, y );
  
axis equal
