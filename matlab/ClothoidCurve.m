classdef ClothoidCurve < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

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
      self.objectHandle = ClothoidCurveMexWrapper( 'new', varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
      %% Destroy the C++ class instance
      ClothoidCurveMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function obj = obj_handle( self )
      obj = self.objectHandle ;
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
      if nargout > 0
        [varargout{1:nargout}] = ClothoidCurveMexWrapper( 'build_G1_D', self.objectHandle, x0, y0, theta0, x1, y1, theta1 );
      else
        ClothoidCurveMexWrapper( 'build_G1', self.objectHandle, x0, y0, theta0, x1, y1, theta1 );
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
      ok = ClothoidCurveMexWrapper( 'build_forward', self.objectHandle, x0, y0, theta0, k0, x1, y1 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = evaluate( self, s )
      % evaluate the curve at curvilinear abscissa `s`
      %
      % Usage:
      %    [x,y] = ref.eval( s )
      %    [x,y,theta,kappa] = ref.eval( s )
      %    
      % On input:
      %    s: curvilinear coordinates where to evaluate the curve
      %       (scalar or vector)
      %
      % On output:
      %    x, y:  coordinates of the curve 
      %    theta: orientation of the curve
      %    kappa: curvature of the curve
      %
      [varargout{1:nargout}] = ClothoidCurveMexWrapper( 'evaluate', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval( self, s )
      [varargout{1:nargout}] = ClothoidCurveMexWrapper( 'eval', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_D( self, s )
      [varargout{1:nargout}] = ClothoidCurveMexWrapper( 'eval_D', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DD( self, s )
      [varargout{1:nargout}] = ClothoidCurveMexWrapper( 'eval_DD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function varargout = eval_DDD( self, s )
      [varargout{1:nargout}] = ClothoidCurveMexWrapper( 'eval_DDD', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [X,Y,S,DST] = closestPoint( self, qx, qy )
      [X,Y,S,DST] = ClothoidCurveMexWrapper( 'closestPoint', self.objectHandle, qx, qy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DST,S] = distance( self, varargin )
      % eval the angle of the circle curve at curvilinear abscissa `s`
      [DST,S] = ClothoidCurveMexWrapper( 'distance', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [X,Y,S,DST] = closestPointBySample( self, qx, qy, ds )
      [X,Y,S,DST] = ClothoidCurveMexWrapper( 'closestPointBySample', self.objectHandle, qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [DST,S] = distanceBySample( self, qx, qy, ds )
      [DST,S] = ClothoidCurveMexWrapper( 'distanceBySample', self.objectHandle, qx, qy, ds );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = xBegin( self )
      res = ClothoidCurveMexWrapper( 'xBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = xEnd( self )
      res = ClothoidCurveMexWrapper( 'xEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yBegin( self )
      res = ClothoidCurveMexWrapper( 'yBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = yEnd( self )
      res = ClothoidCurveMexWrapper( 'yEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaBegin( self )
      res = ClothoidCurveMexWrapper( 'thetaBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = thetaEnd( self )
      res = ClothoidCurveMexWrapper( 'thetaEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = kappaBegin( self )
      res = ClothoidCurveMexWrapper( 'kappaBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = kappaEnd( self )
      res = ClothoidCurveMexWrapper( 'kappaEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = kappa_D( self )
      res = ClothoidCurveMexWrapper( 'kappa_D', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = sMin( self )
      res = ClothoidCurveMexWrapper( 'sMin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = sMax( self )
      res = ClothoidCurveMexWrapper( 'sMax', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = length( self )
      res = ClothoidCurveMexWrapper( 'length', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function info( self )
      fprintf('x0     = %g\n',self.xBegin());
      fprintf('y0     = %g\n',self.yBegin());
      fprintf('theta0 = %g\n',self.thetaBegin());
      fprintf('kappa0 = %g\n',self.kappaBegin());
      fprintf('dk     = %g\n',self.kappa_D());
      fprintf('length = %g\n',self.length());
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function trim( self, smin, smax )
      % trim the clothoid curve at the corresponging curvilinear coordinates
      %
      % Usage:
      %    ref.trim(smin, smax)
      %    
      % On input:
      %    smin:   initial curvilinear coordinate of the curve
      %    smax:   final curvilinear coordinate of the curve
            
      ClothoidCurveMexWrapper( 'trim', self.objectHandle, smin, smax );
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
      ClothoidCurveMexWrapper( 'changeCurvilinearOrigin', self.objectHandle, s0, L );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function rotate( self, angle, cx, cy )
      % rotate the clothoid curve by angle respect to the centre `(cx,cy)`
      %
      % Usage:
      %    ref.rotate(angle, cx, cy)
      %    
      % On input:
      %    angle: the angle of rotation
      %    cx, cy: coordinates of the centre of rotation
      %
      ClothoidCurveMexWrapper( 'rotate', self.objectHandle, angle, cx, cy );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function translate( self, tx, ty )
      % translate the clothoid curve by `(tx,ty)`
      %
      % Usage:
      %    ref.translate(tx, ty)
      %    
      % On input:
      %    tx, ty: horizontal and vertical translation
      %   
      ClothoidCurveMexWrapper( 'translate', self.objectHandle, tx, ty );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeOrigin( self, newX0, newY0 )
      % move the origin of the clothoid to `(newX0, newY0)` 
      %
      % Usage:
      %    ref.changeOrigin(newX0, newY0)
      %    
      % On input:
      %    newX0, newY0: new coordinates of initial point
      %
      ClothoidCurveMexWrapper( 'changeOrigin', self.objectHandle, newX0, newY0 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function scale( self, s )
      % scale clothoid by `sc` factor
      %
      % Usage:
      %    ref.scale(newX0, newY0)
      %    
      % On input:
      %    newX0, newY0: new coordinates of initial point
            
      ClothoidCurveMexWrapper( 'scale', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function reverse( self )
      % reverse the orientation of the clothoid curve 
      % Usage:
      %    ref.reverse()
      %
      ClothoidCurveMexWrapper( 'reverse', self.objectHandle );
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
    function [x0,y0,theta0,k0,dk,L] = getPars( self )
      [x0,y0,theta0,k0,dk,L] = ClothoidCurveMexWrapper( 'getPars', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function BB = bbox( self, max_angle, max_size, varargin )
      % point at infinity
      % Usage:
      %    ref.reverse()
      %
      BB = ClothoidCurveMexWrapper( 'bbox', self.objectHandle, max_angle, max_size, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, npts, varargin )
      % plot: method to plot the clothoid curve
      % Usage:
      %    lineH = ref.plot()
      %    lineH = ref.plot(0.2)
      %    lineH = ref.plot(0.2, 'r', 'LineWidth', 3)
      %    
      % On input:
      %    step:     the distance between consecutive samples used to
      %              draw the curve
      %    varargin: optional arguments passed to the MATLAB plot
      %              function
      %
      % On output:
      %    lineH: the handle to the line object 
      %           (i.e. the output of the MATLAB plot function)
      if nargin < 2
        npts = 1000 ; 
      end
      L     = ClothoidCurveMexWrapper( 'length', self.objectHandle );
      S     = 0:L/npts:L ;
      [X,Y] = ClothoidCurveMexWrapper( 'eval', self.objectHandle, S );
      plot(X,Y, varargin{:});
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotNormal( self, step, len )
      for s=0:step:self.length()
        [x,y,theta,~] = self.evaluate(s) ;
        n = [sin(theta),-cos(theta)] ;
        A = [x,x+len*n(1)];
        B = [y,y+len*n(2)];
        plot(A,B);
      end
    end
  end
end