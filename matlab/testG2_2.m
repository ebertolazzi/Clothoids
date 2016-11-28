clear all ;
close all ;

x0     = -2 ; 
y0     =  3 ;
theta0 = pi/3;
kappa0 = -10.2/1000;

x1     =  3 ;
y1     = 2 ;
theta1 = pi/10; %0*pi/4;
kappa1 = -0.5;

DST = hypot(x1-x0,y1-y0) ;

f = 0.25/DST ;

f0 = 1/4 ;
if abs(kappa0)>8*pi*f
  f0 = 2*pi*f/abs(kappa0) ;
end
f1 = 1/4 ;
if abs(kappa1)>8*pi*f
  f1 = 2*pi*f/abs(kappa1) ;
end
f0
f1
[ S0, S1, SM, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0,x1,y1,theta1,kappa1,f1) ;
S0.L+S1.L+SM.L
flg
draw3curve( S0, S1, SM, true );