clear all ;
close all ;

k0 = -1.484519846508199 ;
k1 = 1.484519846508200 ;

x0     = -2 ; 
y0     =  3 ;
theta0 = pi/2+pi/4;
kappa0 = 10.2/10;

x1     = 3 ;
y1     = 2 ;
theta1 = pi/2; %0*pi/4;
kappa1 = 1.5/4 ;

[ S0, S1, SM, f0, f1, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
for kkk=1:7
  subplot(3,3,kkk) ;
  draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
  title(SM(kkk).opt) ;

  if kkk == 7
    subplot(3,3,8) ;
    draw3angle( S0(kkk), S1(kkk), SM(kkk), true );

    subplot(3,3,9) ;
    draw3curvature( S0(kkk), S1(kkk), SM(kkk), true );
  end
end

%f00 = 0.05 ;
%f11 = 0.05 ;
%[ S0, S1, SM, f0, f1, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f00,x1,y1,theta1,kappa1,f11) ;
%subplot(3,3,7) ;
%draw3curve( S0, S1, SM, true );
