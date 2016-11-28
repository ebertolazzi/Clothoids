close all ;

x0     = -2 ; 
y0     =  3 ;
theta0 = pi/3;
kappa0 = -10.2/10;

x1     =  3 ;
y1     = 2 ;
theta1 = pi/10; %0*pi/4;
kappa1 = -0.5/10;

DST = hypot(x1-x0,y1-y0) ;

close all ;
figure(1);
amin = 0.01 ;
amax = 0.49 ;
astep = (amax-amin)/20 ;
Qmin = 1e100 ;
for f0=amin:astep:amax
  for f1=amin:astep:amax
    fprintf('f0 = %g, f1 = %g\n', f0, f1 ) ;
    [ S0, S1, SM, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0,x1,y1,theta1,kappa1,f1) ;
    if flg
      LEN = S0.L + SM.L + S1.L ;
      if LEN > 5*DST
        fprintf('\nDROPPED\n') ;
        flg = false ;
      end
    end
    if flg
      [tv0,curv0,jerk0] = targetClothoid( S0 ) ;
      [tv1,curv1,jerk1] = targetClothoid( S1 ) ;
      [tvm,curvm,jerkm] = targetClothoid( SM ) ;
      LEN = S0.L + SM.L + S1.L ;
      beta = 2 ;
      %Q = (abs(S0.dkappa)^beta*S0.L+abs(S1.dkappa)^beta*S1.L+abs(SM.dkappa)^beta*SM.L)/LEN ;
      %Q = tv0+tv1+tvm ;
      %Q = curv0+curv1+curvm ;
      Q = (curv0+curv1+curvm)*(tv0+tv1+tvm) ;
      if ( Q < Qmin ) 
        Qmin  = Q ;
        f0min = f0 ;
        f1min = f1 ;
      end
      fprintf('Q = %g\n', Q ) ;
      draw3curve( S0, S1, SM, false );
    else
      fprintf('\nNO OK\n\n') ;
    end
    %Ltot = S0.L+S1.L+2*SM.L;
    %S0, S1, SM, Ltot;
  end
end

x0
y0
theta0
kappa0
f0min
x1
y1
theta1
kappa1
f1min


[ S0, S1, SM, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0min,x1,y1,theta1,kappa1,f1min) ;
draw3curve( S0, S1, SM, true );


f0 = 1/3 ;
if abs(kappa0)>pi/12
  f0 = pi/(4*abs(kappa0)) ;
end
f1 = 1/3 ;
if abs(kappa1)>pi/12
  f1 = pi/(4*abs(kappa1)) ;
end
[ S0, S1, SM, flg ] = buildClothoid3arcG2(x0,y0,theta0,kappa0,f0,x1,y1,theta1,kappa1,f1) ;
draw3curve( S0, S1, SM, true );

