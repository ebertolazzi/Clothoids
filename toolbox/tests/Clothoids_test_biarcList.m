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
  close all;

  %
  % Choosing nodes in parametric curve interpolation
  % E T Y Lee
  %
  % computer-aided design
  % volume 21 number 6 july/august 1989

  PSET = {};

  PSET{1} = [ 0, 0; 26, 24; 28, 24; 54, 0 ];
  PSET{2} = [ 0, 0; 9, 39; 10, 40; 13, 40];
  PSET{3} = [ 0, -0.15; 9.2, -0.15; 10, 0; 10, 0.5831];
  PSET{4} = [ 0, -0.15; 9.2, -0.15; 10, 0; 10.5, 0.3];
  PSET{5} = [ 0, -0.15; 9.2, -0.15; 10, 0; 9.5, 0.3];
  PSET{6} = [ 0, 0; 10, 25; 10, 24; 11, 24.5; 33, 25];
  PSET{7} = [ 0, 41; 10, 41; 31, 41; 40.8, 41; 41, 40.8; ...
              41, 31; 41, 10; 41, 0 ];
  PSET{8}  = [ 0, 41; 40.95, 41; 41, 40.95; 41, 0];
  PSET{9}  = [ 0, 41; 40.95, 41; 41, 40.95; 41, 0; 31, 41; 41, 31];
  PSET{10} = [ 0, 41; 40.95, 41; 41, 40.95; 41, 0; 31, 41; 41, 31; ...
               10, 41; 21, 41; 41, 21; 41, 10];
  PSET{11} = [ 0, 41; 40.95, 41; 41, 40.95; 41, 0; 31, 41; 41, 31; ...
               10, 41; 21, 41; 41, 21; 41, 10; 35, 41; 41, 35];
  PSET{12} = [ 0, 0; 1.34, 5; 5, 8.66; 10, 10; 10.6, 10.4; 10.7, 12; ...
               10.7, 28.6; 10.8, 30.2; 11.4, 30.6; 19.6, 30.6; ...
               20.2, 30.2; 20.3, 28.6; 20.3, 12; 20.4, 10.4; ...
               21, 10; 26, 8.66; 29.66, 5; 31, 0 ];
  PSET{13} = [ 0, 10; 2, 10; 3, 10; 5, 10; 6, 10; 8, 10; 9, 10.5; 11, 15; ...
               12, 50; 14, 60; 15, 85];
  PSET{14} = [ 0, 1; 1, 1.1; 2, 1.1; 3, 1.2; 4, 1.3; 5, 7.2; 6, 3.1; ...
               7, 2.6; 8, 1.9; 9, 1.7; 10, 1.6];
  PSET{15} = [ 2.5, 6.875; 5, 2.23; 10, 0.751; 15, 0.416; 20, 0.283; ...
               25, 0.219; 30, 0.182; 40, 0.143];
  set(0,'DefaultFigureWindowStyle','docked')
  for k=1:15
    figure('Position',[ 1 1 800 800]);
    
    P = PSET{k};
    X = P(:,1);
    Y = P(:,2);

    BL = BiarcList();
    [ theta_guess, theta_min, theta_max, omega, d ] = XY_to_angle( X, Y );

    options = optimoptions('fmincon',...
      'SpecifyObjectiveGradient',false, ...
      'SpecifyConstraintGradient',false, ...
      'CheckGradients', true, ...
      'OptimalityTolerance', 1e-8, ...
      'FunctionTolerance', 1e-8, ...
      'Display','iter' ...
    );

    fun   = @(theta) Target( theta, omega, d );
    theta = fmincon( fun, theta_guess, [], [], [], [], theta_min, theta_max, [], options );
    BL.build_G1( X, Y, theta );
    BL.plot();
    plot( X, Y, 'o', 'LineWidth',2, 'MarkerSize',10,...
          'MarkerEdgeColor','blue', 'MarkerFaceColor',[0.9,0.9,0.9] );
    BL.plotNormal(0.25,0.25);
    title('Biarc List test');
    axis equal;
  end
  set(0,'DefaultFigureWindowStyle','normal')
end
%%%
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