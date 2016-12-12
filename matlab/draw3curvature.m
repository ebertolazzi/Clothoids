function draw3curve(S0,S1,SM,flg)
  % compute points on clothoid
  LEN = S0.L + SM.L + S1.L ;

  X = [ 0, S0.L ]./LEN ;
  K = [ S0.k, S0.k+S0.L*S0.dk ] ;
  plot( X, K, '-r', 'LineWidth', 1.5 ) ;
  hold on

  X = (S0.L+[ 0, SM.L ]) ./ LEN ;
  K = [ SM.k, SM.k+SM.L*SM.dk ] ;
  plot( X, K, '-g', 'LineWidth', 1.5 ) ;

  X =  (S0.L+SM.L+[ 0, S1.L ]) ./ LEN ;
  K = [ S1.k, S1.k+S1.L*S1.dk ] ;
  plot( X, K, '-b', 'LineWidth', 1.5 ) ;

  title('Curvature')
  grid on ;

end
