clear all ;
close all ;

x0     = -1 ;
y0     = 0  ;
theta0 = 2.967059728390360 ;
kappa0 = 0.018315638888734*10 ;

x1     = 1 ;
y1     = 0 ;
theta1 = -3.054326190990077 ;
kappa1 = 0.018315638888734*10 ;

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
