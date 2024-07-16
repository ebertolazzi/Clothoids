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

k_max  = 0.9;
d      = 1/4;
x0     = -d;
y0     = 0;
x3     = d;
y3     = 0;
theta0 = -0*pi/2;
theta3 = linspace(-pi,pi,100);
len    = [];
DlenA  = [];
DlenB  = [];
kind   = [];
subplot(2,2,1);
title('curves');
hold on;
for th=theta3
  DB.build( x0, y0, theta0, x3, y3, th, k_max );
  [l,Da,Db]    = DB.length();
  len(end+1)   = l;
  DlenA(end+1) = Da;
  DlenB(end+1) = Db;
  kind(end+1)  = DB.curve_type();
  DB.plot();
end

subplot(2,2,2);
plot( theta3, len, '-', 'Color', 'blue', 'LineWidth', 2);
title('length vs second angle');

subplot(2,2,3);
hold off
%plot( theta3, DlenA, '-', 'Color', 'blue', 'LineWidth', 2 );
plot( theta3, DlenB, '-', 'Color', 'red', 'LineWidth', 2 );
hold on
if true
  DDlen   = diff(len)./diff(theta3);
  atheta3 = (theta3(1:end-1)+theta3(2:end))/2;
  idx     = find( abs(DDlen) > 100 );
  DDlen(idx) = 0;
  plot( atheta3, DDlen, '-.', 'Color', 'black', 'LineWidth', 2 );
end
ylim([-1.5,1.5]);
legend('analitic','numeric');
title('length derivative');

subplot(2,2,4);
hold off
plot( theta3, kind, '-', 'Color', 'blue');
title('curve type');

