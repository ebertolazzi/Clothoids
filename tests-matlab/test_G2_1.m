addpath('../matlab');

close all ;

x0     = -2 ; 
y0     =  3 ;
theta0 = pi/3;
kappa0 = -0.2;

x1     =  3 ;
y1     = 2 ;
theta1 = -pi/4;
kappa1 = 0.5;

S  = ClothoidList();
ok = S.build3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;

figure(1);
S.plot();
       

