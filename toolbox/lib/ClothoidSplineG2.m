%>
%> Construct a piecewise clothoids \f$ \G(s) \f$ composed by
%> n clothoids that solve the G2 problem
%>
%> **match points**
%>
%> \f[ G(s_k) = \mathbf{p}_k \f]
%> \f[ \lim_{s\to s_k^+} G'(s) = \lim_{s\to s_k^-} G'(s) \f]
%> \f[ \lim_{s\to s_k^+} G''(s) = \lim_{s\to s_k^-} G''(s) \f]
%>
%> **Reference**
%>
%> The solution algorithm is described in
%>
%> - **E.Bertolazzi, M.Frego**, Interpolating clothoid splines with curvature continuity
%>   Mathematica Methods in the Applied Sciences, vol 41, N.4, pp. 1723-1737, 2018
%>
%> \rst
%>
%>   .. image:: ../../images/G2problem3arc.jpg
%>      :width: 80%
%>      :align: center
%>
%> ```
%>
classdef ClothoidSplineG2 < matlab.mixin.Copyable
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    m_objectHandle; % Handle to the underlying C++ class instance
    m_call_delete;
    m_use_PIPAL;
    m_iter_opt;
    m_is_octave;
  end

  methods(Access = protected)
    % make deep copy for copy command
    function obj = copyElement( self )
      obj                = copyElement@matlab.mixin.Copyable(self);
      obj.m_objectHandle = ClothoidSplineG2MexWrapper( 'copy', self.objectHandle );
      obj.m_call_delete  = true;
    end
  end

  methods (Hidden = true)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ceq,jaceq] = nlsys( self, theta )
      ceq   = ClothoidSplineG2MexWrapper( 'constraints', self.m_objectHandle, theta );
      jaceq = ClothoidSplineG2MexWrapper( 'jacobian',    self.m_objectHandle, theta );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [c,ceq,jac,jaceq] = con( self, theta )
      c     = zeros(0,0);
      ceq   = ClothoidSplineG2MexWrapper( 'constraints', self.m_objectHandle, theta );
      jac   = sparse(0,0);
      jaceq = ClothoidSplineG2MexWrapper( 'jacobian', self.m_objectHandle, theta ).'; %'
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [o,g] = obj( self, theta )
      o = ClothoidSplineG2MexWrapper( 'objective', self.m_objectHandle, theta );
      g = ClothoidSplineG2MexWrapper( 'gradient',  self.m_objectHandle, theta );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = ip_constraint( self, theta )
      ceq = ClothoidSplineG2MexWrapper( 'constraints', self.m_objectHandle, theta );
      res = [ceq;-ceq];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = ip_constraint_jacobian( self, theta )
      Jeq = ClothoidSplineG2MexWrapper( 'jacobian', self.m_objectHandle, theta );
      J   = [Jeq;-Jeq];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Call C++ setup for the problem solver.
    %> Check that consecutive points are distinct
    %>
    function build( self, x, y )
      chk = diff(x).^2+diff(y).^2;
      [mi,idx1] = min(chk);
      [ma,idx2] = max(chk);
      if mi == 0 || mi < 1e-10*ma
        error( ...
          [ 'ClothoidSplineG2, build failed: minimum distance between ', ...
            'two consecutive points [min dist=%g, max dist=%g] is 0 or too small!\n', ...
            'index of points at minimum distance is %d\n'], mi, ma, idx1 ...
        );
      end
      ClothoidSplineG2MexWrapper( 'build', self.m_objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = build_internal( self, x, y, varargin )
      if self.m_use_PIPAL
        theta = ClothoidSplineG2MexWrapper( 'pipal', self.m_objectHandle, x, y, varargin{:} );
      else
        %
        % Compute guess angles
        %
        self.build( x, y );
        [ theta_guess, theta_min, theta_max ] = self.guess();
        [~,nc] = self.dims();

        % 'interior-point'
        if self.m_is_octave
          options.TolX   = 1e-10;
          options.TolFun = 1e-20;
        else
          options = optimoptions(...
            'fmincon','Display',self.m_iter_opt, ...
            'CheckGradients',false, ...
            'FiniteDifferenceType','central', ...
            'Algorithm','sqp',...
            'SpecifyConstraintGradient',true,...
            'SpecifyObjectiveGradient',true,...
            'OptimalityTolerance',1e-20,...
            'ConstraintTolerance',1e-10 ...
          );
        end
        obj   = @(theta) self.obj(theta);
        con   = @(theta) self.con(theta);
        theta = fmincon(obj,theta_guess,[],[],[],[],theta_min,theta_max,con,options);
        %options = optimset(varargin{:});
        %[theta,resnorm,~,~,output,~,~] = lsqnonlin( @target, theta, [], [], options );
      end
      %
      % Compute spline parameters
      %
      clots = ClothoidList();
      N     = length(theta);
      clots.reserve(N-1);
      for j=1:N-1
        clots.push_back_G1( x(j), y(j), theta(j), x(j+1), y(j+1), theta(j+1) );
      end
      theta.'
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function clots = build_internal2( self, x, y, varargin )
      %
      % Compute guess angles
      %
      self.build( x, y );
      [ theta_guess, ~, ~ ] = self.guess();
      % 'interior-point'
      if self.m_is_octave
        options.TolX = 1e-20;
      else
        options = optimoptions( ...
          'fsolve', 'Display', self.m_iter_opt, ...
          'CheckGradients', false, ...
          'FiniteDifferenceType', 'central', ...
          'Algorithm', 'levenberg-marquardt',...
          'SpecifyObjectiveGradient', true,...
          'OptimalityTolerance', 1e-20 ...
        );
      end
      obj = @(theta) self.nlsys(theta);
      [theta,~,exitflag,~] = fsolve(obj,theta_guess,options);
      if exitflag <= 0
        error('ClothoidSplineG2, fsolve failed exitflag = %d\n',exitflag);
      end
      %
      % Compute spline parameters
      %
      clots = ClothoidList();
      N     = length(theta);
      clots.reserve(N-1);
      for j=1:N-1
        clots.push_back_G1( x(j), y(j), theta(j), x(j+1), y(j+1), theta(j+1) );
      end
    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = ClothoidSplineG2()
      self.m_objectHandle = ClothoidSplineG2MexWrapper( 'new' );
      self.m_use_PIPAL    = false;
      self.m_iter_opt     = 'iter';
      self.m_is_octave    = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
      if self.m_objectHandle ~= 0
        if self.m_call_delete
          ClothoidSplineG2MexWrapper( 'delete', self.m_objectHandle );
          self.m_objectHandle = 0; % avoid double destruction of object
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function verbose( self, yes )
      if yes
        self.m_iter_opt = 'iter';
      else
        self.m_iter_opt = 'none';
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function use_pipal( self, yesno )
      self.m_use_PIPAL = yesno;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [n,nc] = dims( self )
      [n,nc] = ClothoidSplineG2MexWrapper( 'dims', self.m_objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ theta, theta_min, theta_max ] = guess( self )
      [ theta, theta_min, theta_max ] = ...
        ClothoidSplineG2MexWrapper( 'guess', self.m_objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points (x(i),y(i))
    %> with initial angle theta0 (radiants) and final angle theta1 (radiants)
    %>
    function clots = buildP1( self, x, y, theta0, theta1 )
      ClothoidSplineG2MexWrapper( 'target', self.m_objectHandle, 'P1', theta0, theta1 );
      %clots = self.build_internal2( x, y );
      clots = self.build_internal( x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points (x(i),y(i))
    %> with cyclic condition (tangent and curvature meet)
    %>
    %> - theta(1) = theta(end) mod 2*pi
    %> - kappa(1) = kappa(end)
    %>
    function clots = buildP2( self, x, y )
      S = ClothoidList();
      S.build_G2_cyclic( x, y );
      clots = S;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points (x(i),y(i))
    %> with assignede initial angle and survature
    %> NB: this target is not reccomended (its unstable), used only for debugging
    %>
    function clots = buildP3( self, x, y, theta0, kappa0 )
      ClothoidSplineG2MexWrapper( 'target', self.m_objectHandle, 'P3' );
      %
      % Compute spline parameters
      %
      clots = ClothoidList();
      N     = length(x);
      clots.reserve(N-1);

      c  = ClothoidCurve();
      ok = c.build_forward( x(1), y(1), theta0, kappa0, x(2), y(2) );
      if ok
        clots.push_back( c );
      else
        warning('buildP3 failed');
      end

      for j=3:N
        theta0 = c.thetaEnd();
        kappa0 = c.kappaEnd();
        ok = c.build_forward( x(j-1), y(j-1), theta0, kappa0, x(j), y(j) );
        if ok
          clots.push_back( c );
        else
          warning('buildP3 failed');
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points \f$ (x_i,y_i)\f$
    %> and minimizing derivative of the curvature al the extrema
    %>
    %> minimize \f$ k'(0)^2 + k'(L)^2 \f$
    %>
    function clots = buildP4( self, x, y )
      ClothoidSplineG2MexWrapper( 'target', self.m_objectHandle, 'P4' );
      clots = self.build_internal( x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points \f$ (x_i,y_i)\f$
    %> and minimizing the length of the first and last segment
    %>
    %> minimize \f$ (s_1-s_0)+(s_n - s_{n-1}) \f$
    %>
    function clots = buildP5( self, x, y )
      ClothoidSplineG2MexWrapper( 'target', self.m_objectHandle, 'P5' );
      clots = self.build_internal( x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points \f$ (x_i,y_i)\f$
    %> and minimizing the total length of the spline
    %>
    %> \f$ L = s_n-s_0 \f$
    %>
    function clots = buildP6( self, x, y )
      ClothoidSplineG2MexWrapper( 'target', self.m_objectHandle, 'P6' );
      clots = self.build_internal( x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points \f$ (x_i,y_i)\f$
    %> and minimizing the integral of the square of the curvature:
    %>
    %> minimize: \f$ \displaystyle \int_0^L k(s)^2 \mathrm{d}s \f$
    %>
    function clots = buildP7( self, x, y )
      ClothoidSplineG2MexWrapper( 'target', self.m_objectHandle, 'P7' );
      clots = self.build_internal( x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points \f$ (x_i,y_i)\f$
    %> and minimizing the integral of the square of the curvature derivative:
    %>
    %> minimize:  \f$ \displaystyle \int_0^L k'(s)^2 \mathrm{d}s \f$
    %>
    function clots = buildP8( self, x, y )
      ClothoidSplineG2MexWrapper( 'target', self.m_objectHandle, 'P8' );
      clots = self.build_internal( x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Compute the clothoid spline passing to the points \f$ (x_i,y_i)\f$
    %> and minimizing the integral of the jerk squared
    %>
    function clots = buildP9( self, x, y )
      ClothoidSplineG2MexWrapper( 'target', self.m_objectHandle, 'P9' );
      clots = self.build_internal( x, y );
    end
    %
  end
end
