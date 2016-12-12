function draw3curve(S0,S1,SM,flg)
  % compute points on clothoid

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
    plot( XY(1,:), XY(2,:), '.b', 'LineWidth', 1 ) ;
    hold on

    XY = pointsOnClothoid( SM, 0:SM.L/400:SM.L ) ;
    plot( XY(1,:), XY(2,:), '-k', 'LineWidth', 1 ) ;

    XY = pointsOnClothoid( S1, 0:S1.L/400:S1.L ) ;
    plot( XY(1,:), XY(2,:), '.r', 'LineWidth', 1 ) ;
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

end
