clear all ;
close all ;

x0     = -1 ;
y0     = 0  ;
theta0 = 0  ;
kappa0 = -10 ;

x1     = 1 ;
y1     = 0 ;
kappa1 = -1.5;

for theta1=-pi:pi/50:pi
  [ S, flg ] = buildClothoid(x0,y0,theta0,x1,y1,theta1) ;
  %f0     = 0.1 ; % 0.25 ;
  %f1     = 0.1 ;
  %[ S0, S1, SM, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0,x1,y1,theta1,kappa1,f1) ;
  [ S0, S1, SM, flg, f0, f1 ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
  if flg < 0
    fprintf('provo ricalcolo\n') ;
    [ S0, S1, SM, flg, f0, f1 ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0,x1,y1,theta1,kappa1,f1/10) ;
  end
  if flg < 0
    flg
    f0
    f1
  end
  hold on
  if flg < 0 || S0.L+S1.L+S.L > 10
  else
    draw3curve( S0, S1, SM, false );
  end
end
