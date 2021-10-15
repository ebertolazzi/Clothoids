%=========================================================================%
%                                                                         %
%  Autors: Enrico Bertolazzi                                              %
%          Department of Industrial Engineering                           %
%          University of Trento                                           %
%          enrico.bertolazzi@unitn.it                                     %
%                                                                         %
%=========================================================================%
% Driver test program to check clothoid computation                       %
%=========================================================================%

function test_2_Omega

  X = [2.9265642,2.6734362,2.5109322,1.9078122,1.1859282,1.9249962, ...
     2.8265562,0.00468420000000025,-2.826567,-1.9437558,-1.1859438, ...
     -1.9062558,-2.501565,-2.6734386,-2.9265642,-2.6187522,-1.1406318, ...
     -0.8968758,-1.4562558,-1.9062558,-0.00468780000000013,1.9078122, ...
     1.4468682,0.8968722,1.1406282,2.6187522, 2.9265642 ];
  Y = [-1.707808758,-1.707808758,-2.367185958,-2.582810358,-2.582810358, ...
     -1.167184758,0.915619242,3.178123242,0.915619242,-1.150000758, ...
     -2.582810358,-2.582810358,-2.393750358,-1.707808758,-1.707808758, ...
     -3.178123242,-3.178123242,-2.989063158,-0.915616758,0.925003242, ...
     2.953123242,0.925003242,-0.915616758,-2.989063158,-3.178123242,-3.178123242, -1.707808758 ];

  close all;

  BL = BiarcList();
  [ theta_guess, theta_min, theta_max, omega, d ] = XY_to_angle( X, Y );

  options = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',false, ...
    'SpecifyConstraintGradient',false, ...
    'CheckGradients', true, ...
    'Display','iter' ...
  );

  fun   = @(theta) Target( theta, omega, d );
  theta = fmincon( fun, theta_guess, [], [], [], [], theta_min, theta_max, [], options );

  figure('Position',[ 1 1 800 800]);

  subplot( 1, 2, 1 );
  BL.build_G1( X, Y, theta_guess );
  BL.plot();
  plot( X, Y, 'o', 'LineWidth',2, 'MarkerSize',10,...
        'MarkerEdgeColor','blue', 'MarkerFaceColor',[0.9,0.9,0.9] );
  BL.plotNormal(0.25,0.25);

  title('Biarc List test');
  axis equal;

  subplot( 1, 2, 2 );
  BL.build_G1( X, Y, theta );
  BL.plot();
  plot( X, Y, 'o', 'LineWidth',2, 'MarkerSize',10,...
        'MarkerEdgeColor','blue', 'MarkerFaceColor',[0.9,0.9,0.9] );
  BL.plotNormal(0.25,0.25);
  title('Biarc List test');
  axis equal;

end

function [ res, G ] = Target( theta, omega, d )
  n = length(theta);
  if nargout == 2
    res = 0;
    G   = zeros(n,1);
    for k=1:n-1
      [r,g]    = T( theta(k), theta(k+1), omega(k), d(k) );
      res      = res + r;
      G(k:k+1) = G(k:k+1) + g(:);
    end
  else
    res = 0;
    for k=1:n-1
      r   = T( theta(k), theta(k+1), omega(k), d(k) );
      res = res + r;
    end
  end
end
%
function [ res, G, H ] = T( thetaL, thetaR, omega, d )
  epsilon = 1e-8*d;
  fL      = 4*omega-3*thetaL-thetaR;
  fR      = 4*omega-thetaL-3*thetaR;
  if nargout > 1
    [L, DL] = smooth_sign(fL,epsilon);
    [R, DR] = smooth_sign(fR,epsilon);
    res = L+R;
    G   = DL*[-3,-1]+DR*[-1,-3];
  else
    L   = smooth_sign(fL,epsilon);
    R   = smooth_sign(fR,epsilon);
    res = L+R;
  end
end
%
function [ s, Ds ] = smooth_sign( x, epsilon )
  switch 2
  case 1
    [ s, Ds ] = smooth_sign1( x, epsilon );
  case 2
    [ s, Ds ] = smooth_sign2( x, epsilon );
  case 3
    [ s, Ds ] = smooth_sign3( x, epsilon );
  case 4
    [ s, Ds ] = smooth_sign4( x, epsilon );
  end
end
%
function [ s, Ds ] = smooth_sign1( x, epsilon )
  ax = abs(x)/epsilon;
  e  = exp(-ax);
  s  = epsilon*(2*log1p(e)+ax);
  if nargout > 1
    e1 = e+1;
    if x > 0
      Ds = (1-e)/e1;
    else
      Ds = (e-1)/e1;
    end
  end
end
%
function [ s, Ds ] = smooth_sign2( x, epsilon )
  ax = abs(x);
  if ax > epsilon/2
    s = ax;
    if nargout > 1
      if x > 0
        Ds = 1;
      else
        Ds = -1;
      end
    end
  else
    s = ax*ax/epsilon+epsilon/4;
    if nargout > 1
      Ds = 2*ax/epsilon;
    end
  end
end
%
function [ s, Ds ] = smooth_sign3( x, epsilon )
  s = hypot( 2*epsilon, x );
  if nargout > 1
    Ds = x/s;
  end
end
%
function [ s, Ds ] = smooth_sign4( x, epsilon )
  ax = abs(x);
  if ax > epsilon
    s = ax-epsilon/2;
    if nargout > 1
      if x > 0
        Ds = 1;
      else
        Ds = -1;
      end
    end
  else
    s = ax*ax/(2*epsilon);
    if nargout > 1
      Ds = ax/epsilon;
    end
  end
end
%