classdef ClothoidCurve < CurveBase
  methods
    %> Create a new C++ class instance for the
    %> clothoid arc object
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref = ClothoidCurve()
    %>    ref = ClothoidCurve( x0, y0, theta0, k0, dk, L )
    %>
    %> ```
    %>
    %> **On input:**
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`:   orientation of the clothoid at initial point
    %> - `k0`:       curvature of the clothoid at initial point
    %> - `dk`:       derivative of curvature respect to arclength
    %> - `L`:        length of curve from initial to final point
    %>
    %> **On output:**
    %>
    %> - `ref`: reference handle to the object instance
    %>
    function self = ClothoidCurve( varargin )
      self@CurveBase( 'ClothoidCurveMexWrapper', 'ClothoidCurve' );
      self.objectHandle = ClothoidCurveMexWrapper( 'new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build the clothoid from known parameters
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref.build( x0, y0, theta0, k0, dk, L )
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`:   orientation of the clothoid at initial point
    %> - `k0`:       curvature of the clothoid at initial point
    %> - `dk`:       derivative of curvature respect to arclength
    %> - `L`:        length of curve from initial to final point
    %>
    function build( self, varargin )
      ClothoidCurveMexWrapper( 'build', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Build the interpolating G1 clothoid arc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref.build_G1( x0, y0, theta0, x1, y1, theta1 )
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`:   orientation of the clothoid at initial point
    %> - `x1`, `y1`: coordinate of final point
    %> - `theta1`:   orientation of the clothoid at final point
    %>
    function varargout = build_G1( self, x0, y0, theta0, x1, y1, theta1 )
      if nargout > 1
        [ varargout{1:nargout} ] = ClothoidCurveMexWrapper( ...
          'build_G1_D', self.objectHandle, x0, y0, theta0, x1, y1, theta1 ...
        );
      else
        [ varargout{1:nargout} ] = ClothoidCurveMexWrapper( ...
          'build_G1', self.objectHandle, x0, y0, theta0, x1, y1, theta1 ...
        );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Build the interpolating clothoid arc fixing initial position angle and curvature
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ok = ref.build_forward( x0, y0, theta0, k0, x1, y1 );
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `x0`, `y0`: coordinate of initial point
    %> - `theta0`:   orientation of the clothoid at initial point
    %> - `k0`:       curvature of the clothoid at initial point
    %> - `x1`, `y1`: coordinate of final point
    %>
    %> **On output:**
    %>
    %> - `ok`: true iff the interpolation was successful
    %>
    function ok = build_forward( self, x0, y0, theta0, k0, x1, y1 )
      ok = ClothoidCurveMexWrapper( ...
        'build_forward', self.objectHandle, x0, y0, theta0, k0, x1, y1 ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Find the closest point to a clothoid curve using the
    %> algorithm described in
    %>
    %>  \rst
    %>  | *E. Bertolazzi, M. Frego*, **Point-Clothoid Distance and Projection Computation**,
    %>  | SIAM Journal on Scientific Computing, 2019, 41(5).
    %>  | https://doi.org/10.1137/18M1200439
    %>  ```
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    [ X, Y, S, DST ] = ref.closest_point( qx, qy );
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `qx`, `qx`: coordinate of the point
    %>
    %> **On output:**
    %>
    %> - `X`, `Y` : coordinate of the pojected point
    %> - `S`      : curvilinear coordinate along the clothoid of the projection
    %> - `DST`    : point clothoid distance
    %>
    function [ X, Y, S, DST ] = closest_point( self, qx, qy )
      [ X, Y, S, ~, ~, DST ] =  ClothoidCurveMexWrapper( ...
        'closest_point', self.objectHandle, qx, qy ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function [ X, Y, S, DST ] = closestPoint( self, qx, qy )
      [ X, Y, S, DST ] = self.closest_point( qx, qy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Find the closest point to a clothoid by sampling points
    %>
    %> ```{matlab}
    %>
    %>    [ X, Y, S, DST ] = ref.closest_point_by_sample( qx, qy, ds );
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `qx`, `qx`: coordinate of the point
    %> - `ds`:       sampling distance  coordinate of the point
    %>
    %> **On output:**
    %>
    %> - `X`, `Y` : coordinate of the pojected point
    %> - `S`      : curvilinear coordinate along the clothoid of the projection
    %> - `DST`    : point clothoid distance
    %>
    function [ X, Y, S, DST ] = closest_point_by_sample( self, qx, qy, ds )
      [ X, Y, S, DST ] = ClothoidCurveMexWrapper( ...
        'closest_point_by_sample', self.objectHandle, qx, qy, ds ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function [ X, Y, S, DST ] = closestPointBySample( self, qx, qy, ds )
      [ X, Y, S, DST ] = self.closest_point_by_sample( qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Find the distance between a point and a clothoid by sampling
    %>
    %> ```{matlab}
    %>
    %>    [ DST, S ] = ref.distanceBySample( qx, qy, ds );
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `qx`, `qx`: coordinate of the point
    %> - `ds`:       sampling distance  coordinate of the point
    %>
    %> **On output:**
    %>
    %> - `S`   : curvilinear coordinate along the clothoid of the projection
    %> - `DST` : point clothoid distance
    %>
    function [ DST, S ] = distanceBySample( self, qx, qy, ds )
      [ DST, S ] = ClothoidCurveMexWrapper( ...
        'distanceBySample', self.objectHandle, qx, qy, ds ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Return curvature derivatve of the clothoid
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    dk = ref.dkappa(s0,L);
    %>
    %> ```
    %>
    function res = dkappa( self )
      res = ClothoidCurveMexWrapper( 'dkappa', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> change the origin of the clothoid curve to curviliear corrdinate `s0`
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    ref.change_curvilinear_origin(s0,L);
    %>
    %> ```
    %>
    %> **On input:**
    %>
    %> - `s0`: curvilinear coordinate of the origin of the new curve
    %> - `L`:  nel length of the curve
    %>
    function change_curvilinear_origin( self, s0, L )
      ClothoidCurveMexWrapper( ...
        'change_curvilinear_origin', self.objectHandle, s0, L ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function changeCurvilinearOrigin( self, s0, L )
      self.change_curvilinear_origin( s0, L );
    end
    %% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %function [ s1, s2 ] = intersect( self, C, varargin )
    %  [ s1, s2 ] = ClothoidCurveMexWrapper( ...
    %    'intersect', self.objectHandle, C.obj_handle(), C.is_type(), varargin{:} ...
    %  );
    %end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> point at infinity
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    [xp,yp,xm,ym] = ref.infinity();
    %>
    %> ```
    %>
    %> - `xp`, `yp`: point at infinity (positive arc)
    %> - `xm`, `ym`: point at infinity (negative arc)
    %>
    function [xp,yp,xm,ym] = infinity( self )
      [xp,yp,xm,ym] = ClothoidCurveMexWrapper( 'infinity', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get clothoid parameters
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    [ x0, y0, theta0, k0, dk, L ] = ref.get_pars();
    %>
    %> ```
    %>
    %> - `x0`, `y0`: initial point of the clothoid arc
    %> - `theta0`:   initial angle of the clothoid arc
    %> - `kappa0`:   initial curvature of the clothoid arc
    %> - `dk`:       curvature derivative
    %> - `L`:        length of the clothoid arc
    %>
    function [ x0, y0, theta0, k0, dk, L ] = get_pars( self )
      x0     = self.x_begin();
      y0     = self.y_begin();
      theta0 = self.theta_begin();
      k0     = self.kappa_begin();
      dk     = self.kappa_D(0);
      L      = self.length();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Export clothoid parameters
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    S = ref.export();
    %>
    %> ```
    %>
    %> - `S.x0`, `S.y0`: initial point of the clothoid arc
    %> - `S.x1`, `S.y1`: final point of the clothoid arc
    %> - `S.theta0`:     initial angle of the clothoid arc
    %> - `S.theta1`:     final angle of the clothoid arc
    %> - `S.kappa0`:     initial curvature of the clothoid arc
    %> - `S.kappa1`:     final curvature of the clothoid arc
    %> - `S.dk`:         curvature derivative
    %> - `S.L`:          length of the clothoid arc
    %>
    function S = export( self )
      S.x0     = self.x_begin();
      S.y0     = self.y_begin();
      S.theta0 = self.theta_begin();
      S.kappa0 = self.kappa_begin();
      S.x1     = self.x_end();
      S.y1     = self.y_end();
      S.theta1 = self.theta_end();
      S.kappa1 = self.kappa_end();
      S.dk     = self.kappa_D(0);
      S.L      = self.length();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Get an optimized sampling of curviliear coordinates on the clothoid arc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>    S = ref.optimized_sample(npts,max_angle,offs);
    %>    S = ref.optimized_sample(npts,max_angle,offs,'ISO');
    %>    S = ref.optimized_sample(npts,max_angle,offs,'SAE');
    %>
    %> ```
    %>
    %> **Input:**
    %>
    %> - `npts`:      total number of sampling points
    %> - `max_angle`: initial angle of the clothoid arc
    %> - `offs`:      offset of the curve
    %> - `ISO`:       use ISO orientation of the normal for offset
    %> - `SAE`:       use SAE orientation of the normal for offset
    %>
    %> **Output:**
    %> - `S`: the vector with the sampled curvilinear coordinates
    %>
    function s = optimized_sample( self, npts, max_angle, offs )
      s = ClothoidCurveMexWrapper( ...
        'optimized_sample', self.objectHandle, npts, max_angle, offs ...
      );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the clothoid arc
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot();
    %>   ref.plot( npts );
    %>   ref.plot( npts, 'Color','blue','Linewidth',2);
    %>
    %> ```
    %>
    %> - `npts`: number of sampling points for plotting
    %>
    function plot( self, npts, varargin )
      if nargin < 2
        npts = 1000;
      end
      S = self.optimized_sample( npts, pi/180, 0 );
      [ X, Y ] = ClothoidCurveMexWrapper( 'eval', self.objectHandle, S );
      plot( X, Y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the clothoid arc with offset
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plot_offs( offs );
    %>   ref.plot_offs( offs, npts );
    %>   ref.pplot_offslot( offs, npts, 'Color','blue','Linewidth',2);
    %>
    %> ```
    %>
    %> - `npts`: number of sampling points for plotting
    %> - `offs`:      offset of the curve
    %>
    function plot_offs( self, offs, npts, varargin )
      S = self.optimized_sample( npts, pi/180, offs );
      [ X, Y ] = ClothoidCurveMexWrapper( 'eval', self.objectHandle, S, offs );
      plot( X, Y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the curvature of the clothoid curve
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plotCurvature( npts );
    %>   ref.plotCurvature( npts, 'Color','blue','Linewidth',2);
    %>
    %> ```
    %>
    %> - `npts`: number of sampling points for plotting
    %>
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot_curvature( self, npts, varargin )
      if nargin < 2
        npts = 1000;
      end
      L = ClothoidCurveMexWrapper( 'length', self.objectHandle );
      S = 0:L/npts:L;
      [ ~, ~, ~, kappa ] = ...
        ClothoidCurveMexWrapper( 'evaluate', self.objectHandle, S );
      plot( S, kappa, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotCurvature( self, npts, varargin )
      self.plot_curvature( npts, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the angle of the clothoid curve
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plotAngle( npts );
    %>   ref.plotAngle( npts, 'Color','blue','Linewidth',2);
    %>
    %> ```
    %>
    %> - `npts`: number of sampling points for plotting
    %>
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot_angle( self, npts, varargin )
      if nargin < 2
        npts = 1000;
      end
      L = ClothoidCurveMexWrapper( 'length', self.objectHandle );
      S = 0:L/npts:L;
      [ ~, ~, theta, ~ ] = ...
        ClothoidCurveMexWrapper( 'evaluate', self.objectHandle, S );
      plot( S, theta, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotAngle( self, npts, varargin )
      self.plot_angle( nots, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Plot the normal of the clothoid curve
    %>
    %> **Usage:**
    %>
    %> ```{matlab}
    %>
    %>   ref.plotNormal( step, len );
    %>
    %> ```
    %>
    %> - `step`: number of sampling normals
    %> - `len`:  length of the plotted normal
    %>
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot_normal( self, step, len )
      for s=0:step:self.length()
        [ x, y, theta, ~ ] = self.evaluate(s);
        n = [sin(theta),-cos(theta)];
        A = [x,x+len*n(1)];
        B = [y,y+len*n(2)];
        plot(A,B);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> \deprecated  whill be removed in future version
    %>
    function plotNormal( self, step, len )
      self.plot_normal( step, len );
    end
  end
end
