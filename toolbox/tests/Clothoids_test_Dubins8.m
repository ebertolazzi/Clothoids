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

close all;
clc;
clear all;

addpath('../lib');
addpath('../bin');
addpath('../tests');

DB3a  = Dubins3p();
DB3b  = Dubins3p();

k_max  = 1;

x0     = 7.24;
y0     = 4.75;
theta0 = 0.95;

xM     = 0.73;
yM     = 1.99;

xf     = 5.97;
yf     = 0.67;
thetaf = 0.63;

DB3a.build( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, 'pattern' );
DB3b.build( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, 'ellipse' );

figure();
subplot(1,3,1);
DB3a.plot();
axis equal
grid on
title('pattern');

subplot(1,3,2);
DB3b.plot();
axis equal
grid on
title('ellipse');

subplot(1,3,3);
DB3a.plot();
hold on
DB3b.plot();
axis equal
grid on
title('ellipse');

pa = DB3a.get_pars;
pb = DB3b.get_pars;



fprintf('[ellipse length] %g\n', pb.L123456 );
fprintf('[pattern length] %g\n', pa.L123456 );
fprintf('[ellipse thetam] %g\n', pb.theta3 );
fprintf('[pattern thetam] %g\n', pa.theta3 );
