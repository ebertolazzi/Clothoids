clear all ;
close all ;

x0     = -1 ;
y0     = 0  ;
theta0 = 0  ;
kappa0 = 0.018315638888734*10 ;

x1     = 1 ;
y1     = 0 ;
kappa1 = 0.018315638888734*10 ;
angles = -pi:pi/36:pi ;

for theta0=angles(2:4:end-1)
  theta0
  for theta1=angles(2:4:end-1)
    [ S0, S1, SM, flg, f0, f1 ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
    if flg < 0
      flg
      f0
      f1
    end
    subplot(2,1,1)
    hold on
    draw3curve( S0(7), S1(7), SM(7), false );
    title(SM(7).opt) ;
    subplot(2,1,2)
    hold on
    draw3curve( S0(6), S1(6), SM(6), false );
    title(SM(6).opt) ;
  end
end
