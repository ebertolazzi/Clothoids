clear all ;
close all ;

k0 = -1.484519846508199 ;
k1 = 1.484519846508200 ;

x0     = -2 ; 
y0     =  3 ;
theta0 = pi/2+pi/4;
kappa0 = 10.2/10;

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


[ S0, S1, SM, f0, f1, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
flg
f0
f1
for kkk=1:7
  subplot(3,3,1) ;
  draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
  subplot(3,3,4) ;
  draw3angle( S0(kkk), S1(kkk), SM(kkk), true );
  subplot(3,3,7) ;
  draw3curvature( S0(kkk), S1(kkk), SM(kkk), true );
end

kkk = 1 ;
subplot(3,3,2) ;
draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
title(SM(kkk).opt) ;

kkk = 2 ;
subplot(3,3,3) ;
draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
title(SM(kkk).opt) ;

kkk = 3 ;
subplot(3,3,5) ;
draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
title(SM(kkk).opt) ;

kkk = 4 ;
subplot(3,3,6) ;
draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
title(SM(kkk).opt) ;

kkk = 5 ;
subplot(3,3,8) ;
draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
title(SM(kkk).opt) ;

kkk = 7 ;
subplot(3,3,9) ;
draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
title(SM(kkk).opt) ;

%[ S0, S1, flg ] = buildClothoid2arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
%flg
%S0
%S1

