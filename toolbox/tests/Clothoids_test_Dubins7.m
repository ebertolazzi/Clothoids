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

DB_A = Dubins();
DB_B = Dubins();
DB3  = Dubins3p();

k_max  = 0.6;
d      = 3;
x0     = -d;
y0     = 0;
xM     = 0;
yM     = 0.5*d;
xf     = 2.1*d;
yf     = 0;
theta0 = -pi/2;
thetaf = 0.0;
len    = [];
DlenA  = [];
DlenB  = [];
kind   = [];
epsilon = 1e-4;

DB3.build( x0, y0, theta0, xM, yM, xf, yf, thetaf, k_max, 'pattern' );

thetaGuess = (atan2(yM-y0,xM-x0) + atan2(yf-yM,xf-xM)) / 2;

% non-derivative algorithm to find the optimal thetaM0 (residual = 0)

L  = @(thetaM) len_Dubins(x0, y0, theta0, xM, yM, thetaM, k_max ) + ...
               len_Dubins(xM, yM, thetaM, xf, yf, thetaf, k_max );
DL = @(thetaM) len_Dubins_DR(x0, y0, theta0, xM, yM, thetaM, k_max ) + ...
               len_Dubins_DL(xM, yM, thetaM, xf, yf, thetaf, k_max );

thetaM0 = pattern_search(L);

DB_A.build( x0, y0, theta0,  xM, yM, thetaM0, k_max );
DB_B.build( xM, yM, thetaM0, xf, yf, thetaf,  k_max );

figure();

subplot(2,2,1);

hold on;
DB_A.plot();
DB_B.plot();
plot([x0,xM,xf],[y0,yM,yf],'o','MarkerSize',15,'MarkerFaceColor','red');
axis equal
grid on

subplot(2,2,2);

hold on;
DB3.plot();
plot([x0,xM,xf],[y0,yM,yf],'o','MarkerSize',15,'MarkerFaceColor','red');
axis equal
grid on


npts     = 1000;
thetas   = linspace(thetaGuess-pi,thetaGuess+pi,npts);
LAB      = zeros(1,npts);
D_LAB    = zeros(1,npts);
D_LAB_FD = zeros(1,npts);

for i =1:npts
  th          = thetas(i);
  LAB(i)      = L(th);
  D_LAB(i)    = DL(th);
  D_LAB_FD(i) = (L(th+epsilon)-L(th-epsilon))/(2*epsilon);
end
IDX = find( abs(D_LAB_FD) > 5 );
D_LAB_FD(IDX) = NaN;


subplot(2,2,3);
plot(thetas,LAB,'LineWidth',3);
hold on
plot(thetaGuess,L(thetaGuess),'o','MarkerSize',15,'MarkerFaceColor','magenta');
plot(thetaM0,L(thetaM0),'o','MarkerSize',15,'MarkerFaceColor','blue');
[min_tmp, min_idx] = min(LAB);
plot(thetas(min_idx),min_tmp,'o','MarkerSize',10,'MarkerFaceColor','green');
grid on


subplot(2,2,4);
plot(thetas,D_LAB,'LineWidth',3);
hold on
plot(thetas,D_LAB_FD,'LineWidth',2);
plot(thetaM0,DL(thetaM0),'o','MarkerSize',10,'MarkerFaceColor','green');
grid on
legend({'DL','min','DL_FD'})


%%
function len = len_Dubins( x0, y0, theta0, xf, yf, thetaf, k_max )
  DB = Dubins();
  DB.build( x0, y0, theta0, xf, yf, thetaf, k_max );
  [len,~,~] = DB.length();
end

%%
function DL = len_Dubins_DL( x0, y0, theta0, xf, yf, thetaf, k_max )
  DB = Dubins();
  DB.build( x0, y0, theta0, xf, yf, thetaf, k_max );
  [~,DL,~] = DB.length();
end

%%
function DR = len_Dubins_DR( x0, y0, theta0, xf, yf, thetaf, k_max )
  DB = Dubins();
  DB.build( x0, y0, theta0, xf, yf, thetaf, k_max );
  [~,~,DR] = DB.length();
end

% function to perform a pattern search to find the optimal thetaM0
function thetaM0 = pattern_search( func )
  % Initialize variables
  theta_max       = 2*pi;
  theta_min       = 0;
  theta_candidate = 0;
  numpts          = 16;
  delta_theta     = (theta_max - theta_min) / numpts;
  min_residual    = Inf;
  while delta_theta > 1e-16
    thetas = linspace(theta_min, theta_max, numpts);
    residuals = zeros(1, numpts);
    for i = 1:length(thetas)
      thetaM0      = thetas(i);
      residuals(i) = func(thetaM0);
    end
    [min_residual, min_idx] = min(residuals);
    theta_candidate = thetas(min_idx);
    theta_min       = theta_candidate-delta_theta;
    theta_max       = theta_candidate+delta_theta;
    delta_theta     = (theta_max - theta_min) / numpts;
  end
  thetaM0 = theta_candidate;
end
