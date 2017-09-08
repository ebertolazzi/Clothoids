clear all ;
close all ;

x0     = -1 ;
y0     = 0  ;
theta0 = 2.8*(-0.8) ;
kappa0 = 7.38 ;

x1     = 1 ;
y1     = 0 ;
theta1 = -2.8*(-0.8) ;
kappa1 = 7.38 ;

%x1     = 0.8744165855 ;
%y1     = 0.4704769404 ;
%theta1 = 6.125000000-2*pi ;
%kappa1 = 2.500000000 ;
%x0     = -0.8744165855 ;
%y0     = -0.4704769404 ;
%theta0 = 6.125000000-2*pi ;
%kappa0 = -2.500000000 ;

[ S0, S1, SM, f0, f1, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
for kkk=1:8
  subplot(3,3,kkk) ;
  draw3curve( S0(kkk), S1(kkk), SM(kkk), true );
  %title(SM(kkk).opt) ;

  if kkk == 8
    %subplot(3,3,8) ;
    %draw3angle( S0(kkk), S1(kkk), SM(kkk), true );

    subplot(3,3,9) ;
    draw3curvature( S0(kkk), S1(kkk), SM(kkk), true );
  end
end

%f00 = 0.05 ;
%f11 = 0.05 ;
%[ S0, S1, SM, f0, f1, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f00,x1,y1,theta1,kappa1,f11) ;
%subplot(3,3,7) ;
%draw3curve( S0, S1, SM, true );
