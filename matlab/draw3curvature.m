function draw3curve(S0,S1,SM,flg)

  cm = [1.000000000000000   0.598431372549020   0.200000000000000];

  % compute points on clothoid
  LEN   = S0.L+S1.L+SM.L ;
  
  X = (S0.L+[ 0, SM.L ]) ./ LEN ;
  K = [ SM.k, SM.k+SM.L*SM.dk ] ;
  if flg
    plot( X, K, '--', 'LineWidth', 2, 'Color', cm  ) ;
  else
    plot( X, K, '-r', 'LineWidth', 2 ) ;
  end
  hold on

  X = [ 0, S0.L ]./LEN ;
  K = [ S0.k, S0.k+S0.L*S0.dk ] ;
  plot( X, K, '-b', 'LineWidth', 2 ) ;

  X =  (S0.L+SM.L+[ 0, S1.L ]) ./ LEN ;
  K = [ S1.k, S1.k+S1.L*S1.dk ] ;
  plot( X, K, '-k', 'LineWidth', 2 )  ;

  title('Curvature')
  grid on ;

end
