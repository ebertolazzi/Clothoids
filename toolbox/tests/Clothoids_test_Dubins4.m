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

DB = Dubins();

k_max = 1;

x0     = 0;
y0     = 0;
x3     = 3;
y3     = 0;
theta0 = 0;

kk = 0;
for theta3=linspace(0,pi,16)
  kk = kk+1;
  subplot(4,4,kk);
  DB.build( x0, y0, theta0, x3, y3, theta3, k_max );
  DB.plot();
  axis equal;
  title(sprintf('k=%d angle=%d',kk-1,theta3));
end
