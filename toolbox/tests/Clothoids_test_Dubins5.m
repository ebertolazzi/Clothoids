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

k_max  = 1;
d      = 1/2;
x0     = -d;
y0     = 0;
x3     = d;
y3     = 0;
theta0 = pi/2;
theta3 = linspace(0,pi,400);
len    = [];
Dlen   = [];
kind   = [];
for th=theta3
  DB.build( x0, y0, theta0, x3, y3, th, k_max );
  [l,~,D]   = DB.length();
  len(end+1)  = l;
  Dlen(end+1) = D;
  kind(end+1) = DB.curve_type();
end

subplot(3,1,1);
plot( theta3, len, '-', 'Color', 'blue');

subplot(3,1,2);
hold off
plot( theta3, Dlen, '-', 'Color', 'blue');
if false
  hold on
  DDlen   = diff(Dlen)./diff(len);
  atheta3 = (theta3(1:end-1)+theta3(2:end))/2;
  idx     = find( abs(DDlen) > 100 )
  DDlen(idx) = 0;
  plot( atheta3, DDlen, '-', 'Color', 'red');
end

subplot(3,1,3);
hold off
plot( theta3, kind, '-', 'Color', 'blue');

