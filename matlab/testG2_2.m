clear all ;
close all ;

k0 = -1.484519846508199 ;
k1 = 1.484519846508200 ;

x0     = -2 ; 
y0     =  3 ;
theta0 = pi/2;
kappa0 = -10.2/5;

x1     =  3 ;
y1     = 2 ;
theta1 = pi/2; %0*pi/4;
kappa1 = 1.5;
f0     = 0.067076988604694/2 ; % 0.25 ;
f1     = 0.173268606792847/2 ;

[ S, flg ] = buildClothoid(x0,y0,theta0,x1,y1,theta1) ;
S

%[ S0, S1, SM, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0,x1,y1,theta1,kappa1,f1) ;
%flg
%subplot(2,1,1) ;
%XY = pointsOnClothoid( S, 0:S.L/400:S.L ) ;
%plot( XY(1,:), XY(2,:), '-g', 'LineWidth', 3 ) ;
%hold on
%draw3curve( S0, S1, SM, true );


[ S0, S1, SM, flg, f0, f1 ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
flg
f0
f1
subplot(3,1,1) ;
draw3curve( S0, S1, SM, true );
subplot(3,1,2) ;
draw3angle( S0, S1, SM, true );
subplot(3,1,3) ;
draw3curvature( S0, S1, SM, true );

[ S0, S1, flg ] = buildClothoid2arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
flg
S0
S1

