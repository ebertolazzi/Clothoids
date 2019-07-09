classdef ClothoidCurve < CurveBase
  methods
    function self = ClothoidCurve( varargin )
      %% Create a new C++ class instance for the clothoid arc object
      % Usage:
      %    ref = ClothoidCurve()
      %    ref = ClothoidCurve( x0, y0, theta0, k0, dk, L )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the clothoid at initial point
      %    k0:     curvature of the clothoid at initial point
      %    dk:     derivative of curvature respect to arclength
      %    L:      length of curve from initial to final point
      %
      %  On output:
      %    ref: reference handle to the object instance
      self@CurveBase( 'ClothoidCurveMexWrapper' );
      self.objectHandle = ClothoidCurveMexWrapper( 'new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function str = is_type( ~ )
      str = 'ClothoidCurve';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, varargin )
      %
      % Build the clothoid from known parameters
      %
      % Usage:
      %    ref.build( x0, y0, theta0, k0, dk, L )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the clothoid at initial point
      %    k0:     curvature of the clothoid at initial point
      %    dk:     derivative of curvature respect to arclength
      %    L :     length of curve from initial to final point
      %
      ClothoidCurveMexWrapper( 'build', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = build_G1( self, x0, y0, theta0, x1, y1, theta1 )
      % Build the interpolating G1 clothoid arc
      %
      % Usage:
      %    ref.build_G1( x0, y0, theta0, x1, y1, theta1 )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the clothoid at initial point
      %    x1, y1: coordinate of final point
      %    theta1: orientation of the clothoid at final point
      %
      if nargout > 1
        [ varargout{1:nargout} ] = ...
          ClothoidCurveMexWrapper( 'build_G1_D', ...
                                   self.objectHandle, ...
                                   x0, y0, theta0, x1, y1, theta1 );
      else
        [ varargout{1:nargout} ] = ...
          ClothoidCurveMexWrapper( 'build_G1', ...
                                   self.objectHandle, ...
                                   x0, y0, theta0, x1, y1, theta1 );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = build_forward( self, x0, y0, theta0, k0, x1, y1 )
      % Build the interpolating clothoid arc fixing initial position angle and curvature
      %
      % Usage:
      %    ok = ref.build_forward( x0, y0, theta0, k0, x1, y1 )
      %
      % On input:
      %    x0, y0: coordinate of initial point
      %    theta0: orientation of the clothoid at initial point
      %    k0:     curvature of the clothoid at initial point
      %    x1, y1: coordinate of final point
      %
      % On output:
      %    ok: true iff the interpolation was successful
      %
      ok = ClothoidCurveMexWrapper( 'build_forward', ...
                                    self.objectHandle, ...
                                    x0, y0, theta0, k0, x1, y1 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function aabb_true( self )
      ClothoidCurveMexWrapper( 'aabb_true', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function aabb_false( self )
      ClothoidCurveMexWrapper( 'aabb_false', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ X, Y, S, DST ] = closestPoint( self, qx, qy )
      [ X, Y, S, ~, ~, DST ] = ...
        ClothoidCurveMexWrapper( 'closestPoint', ...
                                 self.objectHandle, qx, qy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ X, Y, S, DST ] = closestPointBySample( self, qx, qy, ds )
      [ X, Y, S, DST ] = ...
        ClothoidCurveMexWrapper( 'closestPointBySample', ...
                                 self.objectHandle, qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ DST, S ] = distanceBySample( self, qx, qy, ds )
      [ DST, S ] = ...
        ClothoidCurveMexWrapper( 'distanceBySample', ...
                                 self.objectHandle, qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = dkappa( self )
      res = ClothoidCurveMexWrapper( 'dkappa', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeCurvilinearOrigin( self, s0, L )
      % change the origin of the clothoid curve to curviliear corrdinate `s0`
      % Usage:
      %    ref.changeOrigin(s0,L)
      %
      % On input:
      %    s0: curvilinear coordinate of the origin of the new curve
      %    L:  nel length of the curve
      %
      ClothoidCurveMexWrapper( 'changeCurvilinearOrigin', ...
                               self.objectHandle, s0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ s1, s2 ] = intersect( self, C, varargin )
      [ s1, s2 ] = ClothoidCurveMexWrapper( 'intersect', ...
                                            self.objectHandle, ...
                                            C.obj_handle(), ...
                                            C.is_type(), ...
                                            varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [xp,yp,xm,ym] = infinity( self )
      % point at infinity
      % Usage:
      %    ref.reverse()
      %
      [xp,yp,xm,ym] = ClothoidCurveMexWrapper( 'infinity', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x0, y0, theta0, k0, dk, L ] = getPars( self )
      x0     = self.xBegin();
      y0     = self.yBegin();
      theta0 = self.thetaBegin();
      k0     = self.kappaBegin();
      dk     = self.kappa_D(0);
      L      = self.length();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [P1,P2,P3] = bbTriangles( self, varargin )
      % bbTriangles( max_angle, max_size, offs );
      [P1,P2,P3] = ClothoidCurveMexWrapper( 'bbTriangles', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function s = optimized_sample( self, npts, max_angle, offs )
      s = ClothoidCurveMexWrapper( 'optimized_sample', ...
                                   self.objectHandle, ...
                                   npts, max_angle, offs );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, npts, varargin )
      if nargin < 2
        npts = 1000;
      end
      S = self.optimized_sample( npts, pi/180, 0 );
      [ X, Y ] = ClothoidCurveMexWrapper( 'eval', self.objectHandle, S );
      plot( X, Y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot_offs( self, offs, npts, varargin )
      S = self.optimized_sample( npts, pi/180, offs );
      [ X, Y ] = ClothoidCurveMexWrapper( 'eval', self.objectHandle, S, offs );
      plot( X, Y, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotCurvature( self, npts, varargin )
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
    function plotTheta( self, npts, varargin )
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
    function plotNormal( self, step, len )
      for s=0:step:self.length()
        [ x, y, theta, ~ ] = self.evaluate(s);
        n = [sin(theta),-cos(theta)];
        A = [x,x+len*n(1)];
        B = [y,y+len*n(2)];
        plot(A,B);
      end
    end
  end
end
