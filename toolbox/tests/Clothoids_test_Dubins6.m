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

DB1 = Dubins();
DB2 = Dubins();
DB3 = Dubins3p();

k_max  = 0.500035;
d      = 3;
x0     = -d;
y0     = 0;
x1     = 1;
y1     = 2;
x2     = d;
y2     = 0;
theta0 = pi/2;
theta1 = linspace(-pi,pi,100);
theta2 = 0;
len    = [];
Dlen   = [];
kind1  = [];
kind2  = [];

%DB3.set_tolerance(1e-20);
DB3.build( x0, y0, theta0, x1, y1, x2, y2, theta2, k_max, 'pattern_trichotomy' );

subplot(2,2,1);
title('curves');
hold on;
for th=theta1
  DB1.build( x0, y0, theta0, x1, y1, th,     k_max );
  DB2.build( x1, y1, th,     x2, y2, theta2, k_max );
  [l1,Da1,Db1] = DB1.length();
  [l2,Da2,Db2] = DB2.length();
  len(end+1)   = l1+l2;
  Dlen(end+1)  = Db1+Da2;
  kind1(end+1) = DB1.curve_type();
  kind2(end+1) = DB2.curve_type();
  DB1.plot();
  DB2.plot();
end


subplot(2,2,2);
plot( theta1, len, '-', 'Color', 'blue', 'LineWidth', 2);
title('length vs second angle');

subplot(2,2,3);
hold off
plot( theta1, Dlen, '-', 'Color', 'red', 'LineWidth', 2 );
hold on
if true
  DDlen   = diff(len)./diff(theta1);
  atheta1 = (theta1(1:end-1)+theta1(2:end))/2;
  idx     = find( abs(DDlen) > 20 );
  DDlen(idx) = NaN;
  plot( atheta1, DDlen, '-.', 'Color', 'black', 'LineWidth', 2 );
end
legend('analitic','numeric');
title('length derivative');

%
% find minimum
%
[~,idx] = min( len );
subplot(2,2,4);
th = theta1(idx);
DB1.build( x0, y0, theta0, x1, y1, th,     k_max );
DB2.build( x1, y1, th,     x2, y2, theta2, k_max );
hold on
DB1.plot();
DB2.plot();
plot( x0, y0, 'o', 'Color', 'blue',  'MarkerSize', 10, 'MarkerFaceColor', 'blue'  );
plot( x1, y1, 'o', 'Color', 'green', 'MarkerSize', 10, 'MarkerFaceColor', 'green' );
plot( x2, y2, 'o', 'Color', 'blue',  'MarkerSize', 10, 'MarkerFaceColor', 'blue'  );
DB3.plot();

axis equal;
title('solution');

DB1.info
DB2.info

DB3.num_evaluation


