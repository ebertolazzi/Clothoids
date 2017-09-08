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
  [ S0, S1, SM, f0, f1, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,x1,y1,theta1,kappa1) ;
  if flg < 0
    flg
    f0
    f1
  else
    subplot(2,1,1) ;
    draw3curve( S0, S1, SM, false );
  end
  [arc1,arc2,ok] = biarc(x0,y0,theta0,x1,y1,theta1);
  if ok
    subplot(2,1,2) ;
    biarc_plot(arc1,arc2,false);
    grid on ;
    axis equal
    hold on
  else
    ok
  end
  %[ S, iter ] = buildClothoid( x0,y0,theta0,x1,y1,theta1 ) ;
  %if iter < 0
  %  iter
  %else
  %  subplot(2,1,2) ;
  %  XY = pointsOnClothoid( S, 0:S.L/400:S.L ) ;
  %  plot( XY(1,:), XY(2,:), '-b', 'LineWidth', 3 ) ;
  %  grid on ;
  %  axis equal
  %  hold on
  %end
end
