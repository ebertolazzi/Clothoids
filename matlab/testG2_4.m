clear all ;
close all ;

k0 = -1.484519846508199 ;
k1 = 1.484519846508200 ;

x0     = -2 ; 
y0     =  3 ;
theta0 = pi/2;
kappa0 = -10.2;

x1     =  3 ;
y1     = 2 ;
theta1 = pi/2; %0*pi/4;
kappa1 = -1.5;
f0     = 0.067076988604694/2 ; % 0.25 ;
f1     = 0.173268606792847/2 ;

hold on
for f0=0.3/50:0.3/50:0.3
  for f1=0.3/50:0.3/50:0.3
    [ S0, S1, SM, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0,x1,y1,theta1,kappa1,f1) ;
    if flg < 0 || S0.L+S1.L+2*SM.L > 20
    else
      flg
      subplot(2,1,1) ;
      draw3curve( S0, S1, SM, false );
      subplot(2,1,2) ;
      draw3angle( S0, S1, SM, false );
    end
  end
end
