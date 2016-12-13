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
  subplot(2,2,1);
  hold on
  draw3curve( S0(7), S1(7), SM(7), false );
  title(SM(7).opt) ;
  
  subplot(2,2,2);
  hold on
  draw3curve( S0(6), S1(6), SM(6), false );
  title(SM(6).opt) ;
  
  subplot(2,2,3);
  hold on
  draw3curve( S0(5), S1(5), SM(5), false );
  title(SM(5).opt) ;
  
  subplot(2,2,4);
  hold on
  draw3curve( S0(4), S1(4), SM(4), false );
  title(SM(4).opt) ;

end
