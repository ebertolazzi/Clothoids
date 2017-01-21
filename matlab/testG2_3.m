clear all ;
close all ;

x0     = -1 ;
y0     = 0  ;
theta0 = 0  ;
kappa0 = -10 ;

x1     = 1 ;
y1     = 0 ;
kappa1 = -1.5;

for theta1=-pi:pi/20:pi
  %[ S, flg ] = buildClothoid(x0,y0,theta0,x1,y1,theta1) ;
  %f0     = 0.1 ; % 0.25 ;
  %f1     = 0.1 ;
  %[ S0, S1, SM, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0,x1,y1,theta1,kappa1,f1) ;
  [ S0, S1, SM, flg, f0, f1 ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
  if flg < 0
    flg
    f0
    f1
  end
  for kkk=1:8
    subplot(2,4,kkk);
    hold on
    draw3curve( S0(kkk), S1(kkk), SM(kkk), false );
    title(SM(kkk).opt) ;
  end
end
