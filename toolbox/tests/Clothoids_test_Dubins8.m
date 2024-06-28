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
DB3c  = Dubins3p();

tol  = 1e-10; % 0.1*m_pi/180.0 };
sang = 2*pi/4;

DB3a.set_tolerance( tol );
DB3a.set_sample_angle( sang );

DB3b.set_tolerance( tol );
DB3b.set_sample_angle( sang );

DB3c.set_sample_points(360);

if false
  k_max  = 0.7497906241084955;

  x0     = -1;
  y0     = 0;
  theta0 = -0.2253391520601129;

  xM     = -0.14345535969001277;
  yM     = -0.14345535969001277;

  xf     = 1;
  yf     = 0;
  thetaf = -0.2253391520601129;
else
  x0     =                             -1;
  y0     =                              0;
  theta0 =         -1.4793802997313432179;
  xM     =        -0.94180274966005206316;
  yM     =        -0.94180274966005206316;
  xf     =                              1;
  yf     =                              0;
  thetaf =         -1.4793802997313432179;
  k_max  =         0.47036903761898174459;
  %thetaM =          4.5902159327450595683; (SAMPLE)
  %len    =          15.846598851367838634; (SAMPLE)
  %thetaM =          4.9447976995266156308; (PATTERN)
  %len    =          16.097284941403650294; (PATTERN)
end

L = @(thetaM) len_Dubins(x0, y0, theta0, xM, yM, thetaM, k_max ) + ...
              len_Dubins(xM, yM, thetaM, xf, yf, thetaf, k_max );

DB3a.build( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, 'pattern_trichotomy' );
DB3b.build( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, 'ellipse' );
DB3c.build( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, 'sample' );

figure();
subplot(2,1,1);
DB3a.plot();
hold on
DB3b.plot();
axis equal
grid on
title('pattern/ellipse');

pa = DB3a.get_pars;
pb = DB3b.get_pars;
pc = DB3c.get_pars;

subplot(2,1,2);

npts   = 1000;
thetas = linspace(0,2*pi,npts);
LAB    = zeros(1,npts);
for i=1:npts
  th     = thetas(i);
  LAB(i) = L(th);
end
plot( thetas, LAB, 'LineWidth', 2 );
hold on;

Lmin = L(pa.theta3);
plot(pa.theta3,Lmin,'o','MarkerSize',12,'MarkerFaceColor','red');
plot([0,2*pi],[Lmin,Lmin],'-','Color','red', 'LineWidth', 2);

samples=DB3a.get_sample_angles( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, tol );
for a=samples
  %plot([a,a],[5,15],'-','LineWidth',2);
  plot(a,L(a),'o','MarkerSize',6,'MarkerFaceColor','yellow');
end

angles=DB3a.get_range_angles( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max );
for a=angles
  plot([a,a],[15,30],'-','LineWidth',2);
  plot(a,L(a),'o','MarkerSize',3,'MarkerFaceColor','red');
end

grid on
title('length');

fprintf('[ellipse length] %.10g\n', pb.L123456 );
fprintf('[pattern length] %.10g\n', pa.L123456 );
fprintf('[sample  length] %.10g\n', pc.L123456 );
fprintf('[ellipse thetam] %.10g\n', pb.theta3  );
fprintf('[pattern thetam] %.10g\n', pa.theta3  );
fprintf('[sample  thetam] %.10g\n', pc.theta3  );

%%
function len = len_Dubins( x0, y0, theta0, xf, yf, thetaf, k_max )
  DB = Dubins();
  DB.build( x0, y0, theta0, xf, yf, thetaf, k_max );
  [len,~,~] = DB.length();
end
