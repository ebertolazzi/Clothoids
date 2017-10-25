function draw3curve(S0,S1,SM,flg)
  LEN   = S0.L+S1.L+SM.L ;
  S     = 0:S0.L/100:S0.L ;
  THETA = S0.theta0+S0.k0.*S+(S0.dk/2).*S.^2 ;
  S     = S./LEN ;
  plot( S, THETA, '-b', 'LineWidth', 3 ) ;
  hold on

  S     = 0:SM.L/100:SM.L ;
  THETA = SM.theta0+SM.k0.*S+(SM.dk/2).*S.^2 ;
  S     = (S0.L+S)./LEN ;
  plot( S, THETA, '-r', 'LineWidth', 3 ) ;

  S     = 0:S1.L/100:S1.L ;
  THETA = S1.theta0+S1.k0.*S+(S1.dk/2).*S.^2 ;
  S     = (S0.L+SM.L+S)./LEN ;
  plot( S, THETA, '-k', 'LineWidth', 3 ) ;

  title('Tangent / Angle')
  grid on ;
end
