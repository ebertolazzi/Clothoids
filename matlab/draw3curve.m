function draw3curve(S0,S1,SM,flg)
  % compute points on clothoid
  
  subplot(2,2,1) ;

  if flg
    XY = pointsOnClothoid( S0, 0:S0.L/400:S0.L ) ;
    plot( XY(1,:), XY(2,:), '-b', 'LineWidth', 3 ) ;
    hold on

    XY = pointsOnClothoid( SM, 0:SM.L/400:SM.L ) ;
    plot( XY(1,:), XY(2,:), '-k', 'LineWidth', 3 ) ;

    XY = pointsOnClothoid( S1, 0:S1.L/400:S1.L ) ;
    plot( XY(1,:), XY(2,:), '-r', 'LineWidth', 3 ) ;
  else
    XY = pointsOnClothoid( S0, 0:S0.L/400:S0.L ) ;
    plot( XY(1,:), XY(2,:), '-r', 'LineWidth', 1 ) ;
    hold on

    XY = pointsOnClothoid( SM, 0:SM.L/400:SM.L ) ;
    plot( XY(1,:), XY(2,:), '-g', 'LineWidth', 1 ) ;

    XY = pointsOnClothoid( S1, 0:S1.L/400:S1.L ) ;
    plot( XY(1,:), XY(2,:), '.b', 'LineWidth', 1 ) ;
  end

  if false
  plot( SM.x, SM.y, '-mo',...
        'LineWidth', 1,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','m',...
        'MarkerSize',5 ) ;

  plot( S0.x, S0.y, '-mo',...
        'LineWidth', 1,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','m',...
        'MarkerSize',5 ) ;

  plot( S1.x, S1.y, '-mo',...
        'LineWidth', 1,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','m',...
        'MarkerSize',5 ) ;

  plot( XY(1,end), XY(2,end), '-mo',...
        'LineWidth', 1,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',5 ) ;
  end
  title('Trajectory')
  grid on ;
  axis equal

  subplot(2,2,2) ;

  if flg
    XY = pointsOnClothoid( S0, 0:S0.L/400:S0.L ) ;
    plot( XY(1,:), XY(2,:), '-b', 'LineWidth', 3 ) ;
    hold on

    XY = pointsOnClothoid( SM, 0:SM.L/400:SM.L ) ;
    plot( XY(1,:), XY(2,:), '-k', 'LineWidth', 3 ) ;

    XY = pointsOnClothoid( S1, 0:S1.L/400:S1.L ) ;
    plot( XY(1,:), XY(2,:), '-r', 'LineWidth', 3 ) ;
    title('Trajectory')
    grid on ;
    axis equal
  end

%  LEN = S0.L + SM.L + S1.L ;

%  X = [ 0, S0.L ]./LEN ;
%  K = [ S0.dkappa, S0.dkappa ] ;
%  plot( X, K, '-r', 'LineWidth', 1.5 ) ;
%  hold on

%  X = (S0.L+[ 0, SM.L ]) ./ LEN ;
%  K = [ SM.dkappa, SM.dkappa ] ;
%  plot( X, K, '-g', 'LineWidth', 1.5 ) ;

%  X =  (S0.L+SM.L+[ 0, S1.L ]) ./ LEN ;
%  K = [ S1.dkappa, S1.dkappa ] ;
%  plot( X, K, '-b', 'LineWidth', 1.5 ) ;

%  title('dCurvature')
%  grid on ;

  subplot(2,2,3) ;

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

  subplot(2,2,4) ;

  S     = 0:S0.L/10:S0.L ;
  THETA = S0.theta+S0.k.*S+(S0.dk/2).*S.^2 ;
  S     = S./LEN ;
  plot( S, THETA, '-r', 'LineWidth', 3 ) ;
  hold on

  S     = 0:SM.L/10:SM.L ;
  THETA = SM.theta+SM.k.*S+(SM.dk/2).*S.^2 ;
  S     = (S0.L+S)./LEN ;
  plot( S, THETA, '-g', 'LineWidth', 3 ) ;

  S     = 0:S1.L/10:S1.L ;
  THETA = S1.theta+S1.k.*S+(S1.dk/2).*S.^2 ;
  S     = (S0.L+SM.L+S)./LEN ;
  plot( S, THETA, '-b', 'LineWidth', 3 ) ;

  title('Tangent / Angle')
  grid on ;
end
